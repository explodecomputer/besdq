# BESDQ — Fast BESD eQTL Query Tool

Fast queries of BESD (Binary Efficient Sequential Data) eQTL summary statistics. Supports direct file access and an optimised SQLite index for repeated queries.

## Why BESDQ?

- **No database server required** — queries run directly against BESD files or a local SQLite index
- **SMR-compatible output** — drop-in replacement for SMR query workflows
- **10–200× faster than direct BESD access** when using the SQLite index for repeated queries
- **Compact storage** — index is 50–70% the size of the original BESD files

See [docs/performance-westra.md](docs/performance-westra.md) for benchmarks.

## Installation

```bash
git clone https://github.com/explodecomputer/besdq.git
cd besdq
pip install -e .
```

## Quickstart

**Query BESD files directly (no setup):**

```bash
besdq --beqtl-summary data/westra_eqtl_hg19 \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output
```

**Query via SQLite index (faster for repeated queries):**

```bash
# One-time index creation
besdq --beqtl-summary data/westra_eqtl_hg19 --index data/westra_eqtl_hg19.db

# Query the index
besdq --besd-index data/westra_eqtl_hg19.db \
  --snp-chrpos 1:100000-2000000 \
  --probe-chrpos 1:1000000-2000000 \
  --out results/output
```

Full usage reference: [docs/usage.md](docs/usage.md)

## Documentation

| Audience | Document |
|---|---|
| Analysts | [docs/usage.md](docs/usage.md) — installation, all query modes, Python API |
| Developers | [docs/developer.md](docs/developer.md) — architecture, file format, design decisions |
| Evaluators | [docs/performance.md](docs/performance.md) — benchmarks |

## Credits

Built on the BESD format from the [SMR software](https://yanglab.westlake.edu.cn/software/smr/#eQTLsummarydata) by [Jian Yang's group](https://yanglab.westlake.edu.cn/).

## License

MIT
