# BESDQ Performance

## Benchmarks

Two benchmark datasets:

| Benchmark | Dataset | Data | Updated |
|---|---|---|---|
| [Westra eQTL](performance-westra.md) | Westra et al. blood eQTL | 506k SNPs × 5,966 probes | Automatically via CI |
| [eQTLGen](performance-eqtlgen.md) | eQTLGen whole-blood cis-eQTL | ~16k traits | Manually (large dataset) |

## Summary (Westra dataset)

| Metric | Direct BESD | SQLite Index |
|---|---|---|
| Query time | 100–200 ms | 10–50 ms |
| Index creation | N/A | ~2 s |
| File size | ~1.2 GB | ~400–600 MB |
| Memory | ~500 MB | ~100 MB |

See [performance-westra.md](performance-westra.md) for detailed results and methodology.
