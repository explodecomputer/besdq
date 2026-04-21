# BESDQ - Fast BESD eQTL Query Tool

Fast queries of BESD (Binary Efficient Sequential Data) eQTL summary statistics files without requiring a database. Supports both SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F formats and provides SMR-compatible command-line interface.

## Installation

### From source

```bash
git clone <repository-url>
cd besdq
pip install -e .
```

This installs the package in development mode with the `besdq` command-line tool.

### Running without installation

You can also run the CLI directly without installing:

```bash
python3 -m besdq.cli --help
```

## Basic Usage

### Creating an SQLite Index Database

For improved performance on repeated queries, create a SQLite index database:

```bash
besdq --beqtl-summary data/westra_eqtl_hg19 --index data/westra_eqtl_hg19.db
```

This creates a database with:
- SQLite tables for SNP and probe metadata with indexed columns for fast range queries
- Per-probe BLOBs storing numpy arrays of associations (snp_indices, betas, SEs)
- Full metadata preservation from the original BESD files

The database can then be used for faster subsequent queries (future implementation).

### Command-line Interface

Query a cis-window using chromosome and position ranges:

```bash
besdq \
  --beqtl-summary data/westra_eqtl_hg19 \
  --query 5e-4 \
  --snp-chr 1 \
  --from-snp-kb 100 \
  --to-snp-kb 2000 \
  --probe-chr 1 \
  --from-probe-kb 1000 \
  --to-probe-kb 2000 \
  --out results/westra_1
```

#### Coordinate Format Options

You can specify positions in multiple ways:

**Kilobase format (default):**
```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chr 1 --from-snp-kb 100 --to-snp-kb 2000 \
  --probe-chr 1 --from-probe-kb 1000 --to-probe-kb 2000 \
  --out results/output
```

**Base pair format:**
```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chr 1 --from-snp-bp 100000 --to-snp-bp 2000000 \
  --probe-chr 1 --from-probe-bp 1000000 --to-probe-bp 2000000 \
  --out results/output
```

**Chr:pos format (range or single position):**
```bash
# Range
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output

# Single position
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:1191870 \
  --probe-chrpos 1:1140818 \
  --out results/output
```

#### P-value Filtering

Filter results by p-value threshold (default: 0.05):

```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --query 1e-4 \
  --out results/output
```

### Python API

Use the query engine directly in Python:

```python
from besdq import BESDQueryEngine

# Initialize
engine = BESDQueryEngine('data/westra_eqtl_hg19')

# Query associations
associations = engine.query_cis_window(
    snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
    probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
)

# Results include SNP/probe metadata and statistics
for assoc in associations:
    print(f"{assoc['snp_id']} - {assoc['probe_id']}: "
          f"beta={assoc['beta']:.4f}, p={assoc['pval']:.2e}")
```

## Output Format

Results are written in tab-separated format compatible with SMR:

```
SNP    SNP_Chr  SNP_bp   A1  A2   Probe           Probe_Chr  Probe_bp  Gene     Beta       SE         P_value
rs123  1        1191870  T   C    ILMN_2349633    1          1140818   TNFRSF18 -0.436080  0.040022   1.23e-25
...
```

## Unit Testing

Run the test suite:

```bash
python3 -m unittest tests.test_queries -v
```

Run a specific test:

```bash
python3 -m unittest tests.test_queries.TestBESDQueryEngine.test_single_position_query -v
```

### Test Coverage

The test suite (`tests/test_queries.py`) includes:
- Data loading verification (SNP/probe counts, format detection)
- Single position queries
- Range queries
- P-value calculation accuracy
- Beta and SE value storage
- SNP and probe metadata indexing
- Chromosome and position filtering
- Empty query handling

All tests use the `westra_eqtl_hg19` dataset as reference data.

## File Format

BESDQ expects three files with a common prefix:

- `.besd` - Binary BESD file with association statistics
- `.esi` - SNP index file (chr, rsid, genetic_distance, bp, allele1, allele2, frequency)
- `.epi` - Probe index file (chr, probe_id, genetic_distance, probe_bp, gene, orientation)

## Requirements

- Python 3.9+
- numpy (for efficient array storage in database indices)
- No external dependencies for core query functionality when not using database indices

## License

MIT
