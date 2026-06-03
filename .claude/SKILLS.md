# BESDQ Project Skills Reference

The canonical developer reference is [docs/developer.md](../docs/developer.md).
The domain glossary is [CONTEXT.md](../CONTEXT.md).
Design decision rationale is in [docs/adr/](../docs/adr/).

## Quick orientation for AI assistants

- **Two query engines**: `BESDQueryEngine` (direct BESD) and `BESDQueryIndex` (SQLite). Both share the same query interface.
- **Statistics encoding**: z-scores stored as float16. VectorN mode (GWAS-SSF) stores SE per SNP; ScalarN mode (legacy BESD) stores n per trait and reconstructs SE.
- **SNP key convention**: `chr:pos:A1:A2` where A1 is alphabetically first. Alleles and beta are normalised at import time.
- **Two-stage build**: Stage 1 downloads + filters per-trait intermediates (resumable); Stage 2 collates them into the SQLite index (always clean rebuild).
- **Output scale**: SD units by default (`beta / sd_y`); `--original-scale` returns original units.
- **`probe`** is the historical BESD synonym for **`trait`** — used interchangeably in internal code.

## Package layout

```
besdq/
├── besd_reader.py     # BESDQueryEngine, BESDReader, IndexReader, norm_cdf
├── sqlite_query.py    # BESDQueryIndex, get_probe_snps
├── builder.py         # BESDIndexBuilder
├── queries.py         # Unified re-export
└── cli.py             # CLI entry point
```

## Common patterns

```python
# Direct query
from besdq import BESDQueryEngine
engine = BESDQueryEngine('data/westra')
assocs = engine.query_cis_window(snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
                                  probe_chr='1', probe_start_kb=1000, probe_end_kb=2000)

# SQLite query
from besdq import BESDQueryIndex
with BESDQueryIndex('data/westra.db') as idx:
    assocs = idx.query_by_snp_id('rs3818646')

# Build index
from besdq import BESDIndexBuilder
BESDIndexBuilder('data/westra.db').build('data/westra', force=True)
```

## Tests

```bash
python3 -m unittest tests.test_queries -v
```

Test data: Westra eQTL (506,049 SNPs, 5,966 probes) in `data/`.
