# BESDQ Developer Reference

Architecture, file format specification, and implementation details for contributors and maintainers. For the *why* behind major decisions, see [docs/adr/](adr/).

---

## Package Structure

```
besdq/
├── __init__.py              # Public API exports
├── besd_reader.py           # BESDQueryEngine — direct BESD file access
│   ├── BESDQueryEngine      # Main query interface
│   ├── BESDReader           # Binary BESD file parser
│   ├── IndexReader          # .esi / .epi file reader
│   └── norm_cdf()           # P-value via Abramowitz-Stegun approximation
├── sqlite_query.py          # BESDQueryIndex — SQLite queries
│   ├── BESDQueryIndex       # Query interface (context manager)
│   └── get_probe_snps()     # Deserialize numpy BLOBs
├── builder.py               # BESDIndexBuilder — database creation
│   ├── build()              # Main build workflow
│   ├── _create_schema()     # SQLite schema setup
│   └── _write_probe_data()  # Serialize numpy arrays to BLOBs
├── queries.py               # Unified interface (re-exports both engines)
└── cli.py                   # Command-line interface
```

---

## Architecture

### Hybrid Storage Design

- **Metadata & indices** (`esi`, `epi` tables): SQLite with indexed `(chr, bp)` columns for O(log n) range queries over millions of SNPs/probes.
- **Statistics** (`probe_data` table): per-probe BLOBs storing binary-serialised numpy arrays. Preserves sequential read locality from the original sparse BESD format.

### SQLite Schema

```sql
CREATE TABLE esi (
  row_idx INTEGER PRIMARY KEY,
  chr TEXT, snp_id TEXT, genetic_dist REAL, bp INTEGER,
  a1 TEXT, a2 TEXT, freq REAL
);
CREATE INDEX idx_esi_chr_bp ON esi(chr, bp);
CREATE INDEX idx_esi_snp_id ON esi(snp_id);

CREATE TABLE epi (
  row_idx INTEGER PRIMARY KEY,
  chr TEXT, probe_id TEXT, genetic_dist REAL, probe_bp INTEGER,
  gene TEXT, orientation TEXT
);
CREATE INDEX idx_epi_chr_bp ON epi(chr, probe_bp);
CREATE INDEX idx_epi_probe_id ON epi(probe_id);

CREATE TABLE probe_data (
  probe_idx INTEGER PRIMARY KEY,
  snp_count INTEGER,
  snp_indices BLOB,  -- int32 numpy array
  zscores BLOB,      -- float16 numpy array
  se_vector BLOB,    -- float16 numpy array (VectorN only; NULL in ScalarN)
  n_scalar INTEGER   -- ScalarN only (NULL in VectorN)
);

CREATE TABLE besd_meta (
  key TEXT PRIMARY KEY,
  value TEXT
);
```

See [ADR-0001](adr/0001-zscore-statistics-encoding.md) for why z-scores are stored as float16 rather than raw betas, and [ADR-0002](adr/0002-vector-mode-stores-se-not-n.md) for the VectorN vs ScalarN split.

---

## BESD File Format (Sparse)

The binary format read by `BESDReader`. Source: [SMR source code](https://yanglab.westlake.edu.cn/software/smr/#Sourcecode) `src/SMR_data_p1.cpp`.

### Magic Numbers

| Value | Format |
|---|---|
| `0x40400000` | `SPARSE_FILE_TYPE_3F` |
| `3` | `SPARSE_FILE_TYPE_3` |

Dense format is not supported.

### File Structure

```
[uint32  magic]
[15 × uint32  reserved]   -- type 3 only, absent in type 3F
[uint64  val_num]          -- total number of stored values
[uint64  cols[2 × n_probes + 1]]  -- per-probe boundary array
[uint32  rowid[val_num]]   -- SNP row indices
[float32 val[val_num]]     -- interleaved betas and SEs
```

### Layout Example

| | probe0 beta | probe0 SE | probe1 beta | probe1 SE |
|---|---|---|---|---|
| snp0 | 0.2 | 0.05 | | |
| snp1 | | | 0.3 | 0.07 |
| snp2 | -0.1 | 0.04 | | |

Serialised left-to-right, betas before SEs per probe:

| index | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| value | 0.2 | -0.1 | 0.05 | 0.04 | 0.3 | 0.07 |
| probe | 0 | 0 | 0 | 0 | 1 | 1 |
| field | beta | beta | SE | SE | beta | SE |
| snp | 0 | 2 | 0 | 2 | 1 | 1 |

- `val_num = 6`
- `cols[5] = [0, 2, 4, 5, 6]` — beta range for probe 0 is `[0,2)`, SE range is `[2,4)`, etc.
- `rowid[6] = [0, 2, 0, 2, 1, 1]`
- `val[6] = [0.2, -0.1, 0.05, 0.04, 0.3, 0.07]`

### Reconstructing an Association

For probe `p`, association `i`:

1. `beta_start = cols[2p]`
2. `se_start = cols[2p+1]`
3. `count = se_start - beta_start`
4. `snp_idx = rowid[beta_start + i]`
5. `beta = val[beta_start + i]`
6. `se = val[se_start + i]`

`val_num` is the global payload length for reading `rowid` and `val`; per-probe counts come from `cols` boundaries.

---

## Statistics Encoding

Two storage modes; selected at build time based on data source:

**VectorN** (GWAS-SSF and any source providing SE directly):
- Stores `zscores` (float16) + `se_vector` (float16, original-unit SE) per SNP-trait pair.
- At query time: `se_out = se_stored / sd_y`, `beta_out = z × se_out`.

**ScalarN** (legacy BESD imports):
- Stores `zscores` (float16) per pair + `n_scalar INTEGER` per trait.
- SE reconstructed from n, AF, and trait_var at query time.

See [ADR-0001](adr/0001-zscore-statistics-encoding.md) and [ADR-0002](adr/0002-vector-mode-stores-se-not-n.md).

### Output Scale

Default output is **SD units**: `beta_sd = beta_orig / sd_y`, `se_sd = se_orig / sd_y` where `sd_y = sqrt(trait_var)`. Z-score and p-value are invariant. Pass `--original-scale` to return original-unit values. See [ADR-0004](adr/0004-vector-mode-and-sd-unit-output.md).

---

## P-value Calculation

```python
pval = 2 * (1 - norm_cdf(abs(beta / se)))
```

Uses the Abramowitz and Stegun polynomial approximation (~0.00012 accuracy). Edge cases:
- SE = 0: returns p = 0 (skips calculation).
- Extreme z-scores: p underflows to 0 — expected for very significant associations.

---

## Numpy Array Serialisation

- Serialise: `array.tobytes()`
- Deserialise: `numpy.frombuffer(blob, dtype=...)`
- Stored dtypes: `int32` for SNP indices, `float16` for z-scores and SEs.
- Convert `numpy.int32` to Python `int` before SQLite parameterised queries.

---

## SNP Identity Convention

The canonical SNP key is `chr:pos:A1:A2` where **A1 is the alphabetically first allele**. At import, if the source effect allele is not alphabetically first, alleles are swapped and beta/z-score is negated. ESI deduplication uses this key. `snp_id` (rsid) is a lookup field, not the primary key; variants without an rsid are stored with `snp_id = NULL`.

---

## Significance Mask (GWAS-SSF Import)

See [ADR-0003](adr/0003-gwas-ssf-significance-mask.md) for rationale.

| Tier | Condition | Stored |
|---|---|---|
| Cis | Within cis-radius of `trait_chr`/`trait_bp` | All variants |
| Significant trans | p < 5×10⁻⁸, outside cis | ±sig-radius window around each LD-independent peak |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ | Dropped |

Trans peak identification uses plink2 LD clumping (r²=0.01, window=10,000 kb by default).

---

## Two-Stage Build

See [ADR-0006](adr/0006-two-stage-build-with-persistent-intermediates.md) for design rationale.

- **Stage 1** (`besdq-build-stage1`): per-trait download → filter → write intermediate pair (`GCST*.tsv.gz` + `GCST*.yaml`) → delete download. Resumable: completed traits are skipped; failed traits (`GCST*.failed.yaml`) are retried.
- **Stage 2** (`besdq-build-stage2`): consume all completed intermediates → build SQLite index. Always a clean rebuild.

---

## Query Modes

### `BESDQueryEngine` (direct BESD)

- Loads `.esi`, `.epi`, and `.besd` arrays into memory on init.
- Query time: ~100–200 ms.
- Use for one-off queries; no setup required.

### `BESDQueryIndex` (SQLite)

- Metadata in indexed SQLite tables; statistics in per-probe BLOBs.
- Query time: ~10–50 ms.
- Context manager handles connection lifecycle.
- Index creation: ~2 s (one-time).

### Supported Query Types

1. Cis-window: SNP region + probe region simultaneously
2. SNP ID: all associations for a specific SNP
3. Probe/trait ID: all associations for a specific probe
4. Gene: all associations for a specific gene

---

## Testing

26 unit tests covering all query types, coordinate formats, p-value accuracy, and consistency between BESD and SQLite modes. Test data: Westra eQTL (506,049 SNPs, 5,966 probes).

```bash
python3 -m unittest tests.test_queries -v
```

---

## Known Limitations

- Result ordering differs between BESD and SQLite modes (same data, different order).
- Dense BESD format is not supported.
- P-values can underflow to 0 for extremely significant associations (z > ~8.3).
