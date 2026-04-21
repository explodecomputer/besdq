**Project: `besdq` — a Python library for querying BESD eQTL summary data**

**Architecture decision:** Hybrid SQLite store. Metadata and indices (ESI, EPI) are stored as normal SQL tables with indexed columns for fast range and ID queries. Summary statistics are stored as per-probe BLOBs (binary-serialised numpy arrays) to preserve sequential read locality from the original sparse BESD format.

**Schema:**
- `esi(row_idx PK, chr, snp_id, genetic_dist, bp, a1, a2, freq)` — indexed on `(chr, bp)` and `snp_id`
- `epi(row_idx PK, chr, probe_id, genetic_dist, probe_bp, gene, orientation)` — indexed on `(chr, probe_bp)` and `probe_id`
- `probe_data(probe_idx PK, snp_count, snp_indices BLOB, betas BLOB, ses BLOB)` — int32/float32 numpy arrays serialised with `tobytes()`
- `besd_meta(key PK, value)` — stores format, n_snps, n_probes, source file paths

**Query modes to support:**
1. Variant range query: `chr:start-end` → SNP row indices via SQL → filter probe BLOBs in numpy
2. Probe range/ID query: probe region or ID → probe row indices via SQL → deserialise BLOBs directly
3. Cis-window query: SNP region + probe region simultaneously — get SNP indices from SQL, get probe indices from SQL, intersect in numpy

**BESD format notes (sparse):**
- Magic number in first 4 bytes distinguishes dense (`0x40000000`) from sparse (`0x3f800000`)
- Sparse header: `[magic][reserved][N_probes][N_snps]` as uint32s
- Per-probe offset table follows: `N_probes × 2 × uint64` (start, end byte offsets into file)
- Each probe block: array of `uint32` SNP indices, then `float32` betas, then `float32` SEs
- Row order in `.esi`/`.epi` must be preserved exactly — these indices are baked into the BESD offsets

**Key implementation files to create:**
- `besd.py` — raw BESD binary parser (reads header, offset table, deserialises probe blocks)
- `builder.py` — ingests `.esi`, `.epi`, `.besd` and populates the SQLite DB
- `query.py` — the public query API against the SQLite DB
- `cli.py` — optional command-line interface wrapping the query API

