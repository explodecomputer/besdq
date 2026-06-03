# BESDQ Usage Guide

## Installation

### From GitHub (recommended)

```bash
git clone https://github.com/explodecomputer/besdq.git
cd besdq
pip install -e .
```

### From PyPI (when available)

```bash
pip install besdq
```

### Without installation

```bash
python3 -m besdq.cli --help
```

## Dependencies

**Required:**
- Python 3.9+
- numpy >= 1.20
- PyYAML >= 5.1

**System requirements:**
- ~1 GB free disk space per BESD dataset
- ~50–70% additional space for the SQLite index (optional)

**For LD clumping during trans-eQTL import:**
- [plink2](https://www.cog-genomics.org/plink/2.0/) on `PATH` (`conda install -c bioconda plink2`)

---

## Two Query Modes

| Mode | Flag | When to use |
|---|---|---|
| Direct BESD | `--beqtl-summary` | One-off queries, exploratory analysis |
| SQLite index | `--besd-index` | Repeated queries, large-scale analysis |

---

## Querying BESD Files Directly

No setup required — works immediately on `.besd` / `.esi` / `.epi` file triples.

```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output
```

---

## SQLite Index

### Creating the index (one-time)

```bash
besdq --beqtl-summary data/westra_eqtl_hg19 --index data/westra_eqtl_hg19.db
```

### Querying the index

```bash
besdq --besd-index data/westra_eqtl_hg19.db \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output
```

---

## Coordinate Formats

All query modes accept the same coordinate options:

**Chr:pos (recommended):**
```bash
# Range
--snp-chrpos 1:100000-2000000
# Single position
--snp-chrpos 1:1191870
```

**Kilobase:**
```bash
--snp-chr 1 --from-snp-kb 100 --to-snp-kb 2000
```

**Base pair:**
```bash
--snp-chr 1 --from-snp-bp 100000 --to-snp-bp 2000000
```

Validation rules:
- For `chr:start-end`, `start` ≤ `end`.
- `--snp`, `--probe`, and `--gene` are mutually exclusive.
- Identifier flags (`--snp`, `--probe`, `--gene`) cannot be combined with region flags.

---

## P-value Filtering

```bash
besdq --besd-index data/westra_eqtl_hg19.db \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --query 1e-4 \
  --out results/output
```

---

## Output Format

Tab-separated, SMR-compatible:

```
SNP    SNP_Chr  SNP_bp   A1  A2  Probe         Probe_Chr  Probe_bp  Gene      Beta       SE         P_value
rs123  1        1191870  T   C   ILMN_2349633  1          1140818   TNFRSF18  -0.436080  0.040022   1.23e-25
```

---

## Building an Index from GWAS-SSF Files (EBI GWAS Catalog)

If your data are in [GWAS-SSF format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format), use `import-gwas-ssf` to build a queryable index directly.

### 1. Prepare a trait annotation file

Create a tab-separated file with one row per summary-statistics file:

```
file_path	trait_id	trait_name	trait_chr	trait_bp	gene	context
data/GCST90275731.h.tsv.gz	GCST90275731	IL10 expression PBMC	1	206774541	IL10	PBMC_baseline
```

Required columns: `file_path`, `trait_id`, `trait_name`.
Optional columns:

| Column | Default | Description |
|---|---|---|
| `trait_chr` | — | Chromosome of trait locus (enables cis storage) |
| `trait_bp` | — | Base-pair position of trait locus |
| `sample_size` | from YAML | Sample size N |
| `trait_var` | 1.0 | Phenotype variance |
| `gene` | — | Gene symbol |
| `context` | — | Experimental context (tissue, treatment, etc.) |

If `sample_size` is omitted, it is read from the companion `*-meta.yaml` file provided by EBI.

### 2. Run the import

```bash
import-gwas-ssf \
  --trait-annotation data/traits.tsv \
  --ld-reference data/ldref/EUR \
  --workers 3 \
  --output data/study.db
```

**Key options:**

| Option | Default | Description |
|---|---|---|
| `--trait-annotation FILE` | — | Trait annotation TSV (required) |
| `--ld-reference PREFIX` | — | plink2 LD reference prefix (required) |
| `--output FILE` | `<tsv_stem>.db` | Output database |
| `--workers N` | 1 | Parallel workers |
| `--cis-radius BP` | 1,000,000 | Cis window radius |
| `--sig-threshold P` | 5e-8 | Genome-wide significance threshold |
| `--sug-threshold P` | 1e-4 | Suggestive threshold |
| `--sig-radius BP` | 500,000 | Trans peak window radius |

**Significance filtering applied during import:**

| Tier | Condition | Stored |
|---|---|---|
| Cis | SNP within ±1 Mb of `trait_chr`/`trait_bp` | All variants |
| Significant trans | p < 5×10⁻⁸, different chromosome | ±500 kb window around each independent peak |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ | Dropped |

### 3. Query

```bash
besdq --db data/study.db --info
besdq --db data/study.db --snp rs12238997 --out results/out
besdq --db data/study.db --probe GCST90275731 --out results/out
```

---

## Large-Scale Builds: Two-Stage Pipeline

For publications with hundreds or thousands of traits, use the two-stage pipeline. Stage 1 is resumable per-trait; peak disk usage is bounded.

### 1. Discover traits for a publication

```bash
besdq-discover-study 38714679 \
  --parse-gene \
  --gene-annotation resources/ensembl_gene_tss_hg38.tsv \
  --output data/38714679/discovery.tsv
```

Queries the EBI GWAS Catalog API for all GCST accessions under the given PMID. With `--parse-gene`, maps gene symbols to chromosomal positions via a local Ensembl TSV (`gene_name`, `chr`, `bp` columns required).

### 2. Stage 1 — download, filter, write intermediates

```bash
besdq-build-stage1 data/38714679/discovery.tsv \
  --workdir data/38714679 \
  --parallel-downloads 4
```

For each trait: streams the GWAS-SSF file, applies the significance filter, writes a `GCST*.tsv.gz` + `GCST*.yaml` pair, deletes the download. Local paths are never deleted.

**Resumability:** re-running skips completed traits, retries failed ones (`GCST*.failed.yaml`), ignores partial writes.

| Option | Default | Description |
|---|---|---|
| `--workdir DIR` | — | Working directory for intermediates (required) |
| `--parallel-downloads N` | 1 | Concurrent download workers |
| `--cis-radius BP` | 1,000,000 | Cis window radius |
| `--ld-reference PREFIX` | — | plink2 LD reference for trans clumping |

### 3. Stage 2 — collate into queryable index

```bash
besdq-build-stage2 data/38714679/38714679 study_38714679.db
```

Always a clean rebuild. Query with `besdq --db study_38714679.db`.

---

## Python API

```python
from besdq import BESDQueryEngine

engine = BESDQueryEngine('data/westra_eqtl_hg19')
assocs = engine.query_cis_window(
    snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
    probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
)
for a in assocs:
    print(f"{a['snp_id']} — {a['trait_id']}: beta={a['beta']:.4f}, p={a['pval']:.2e}")
```

```python
from besdq import BESDQueryIndex

with BESDQueryIndex('data/westra_eqtl_hg19.db') as index:
    assocs = index.query_by_probe_id('ILMN_2349633')
    assocs = index.query_by_snp_id('rs3818646')
    assocs = index.query_cis_window(
        snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
        probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
    )
```

---

## Running Tests

```bash
python3 -m unittest tests.test_queries -v
```
