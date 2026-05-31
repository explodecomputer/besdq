# ADR 0006 — Two-stage build with persistent intermediate files

**Status**: Accepted

## Context

Building a BESDQ index from thousands of EBI GWAS-SSF files requires downloading files of ~500MB–2GB each. The current single-pass pipeline holds filtered rows in memory and writes the index in one run. At scale this is impractical: a crash mid-run loses all progress, and keeping thousands of downloads on disk simultaneously is prohibitive.

## Decision

Split the build into two user-facing stages with a persistent on-disk intermediate between them.

**Stage 1** (`build-stage1`): download → filter → write intermediate → delete download. One intermediate file pair per trait (`GCST<id>.tsv.gz` + `GCST<id>.yaml`), written atomically via `.tmp` rename. Intermediates are grouped under `<workdir>/<pmid>/`. Failed traits write `GCST<id>.failed.yaml` and are retried on re-run; completed traits are skipped. Bounded-concurrency download via `--parallel-downloads N` (each worker owns the full download → filter → delete lifecycle).

**Stage 2** (`build-stage2`): read all completed intermediates in a PMID directory, collate into SQLite index. Always a clean rebuild — no resumability.

Intermediate format is gzipped TSV + YAML sidecar (human-readable, no new binary dependencies).

## Alternatives considered

**Single-pass with checkpointing**: checkpoint the in-memory state to disk periodically. Rejected — checkpointing in-memory filter state is complex and the intermediate is not independently inspectable.

**Parquet intermediates**: more compact and faster to read. Rejected — requires `pyarrow`, not human-readable, and the volume of intermediate data is modest (already filtered to cis/trans only).

**Parallel download with unbounded concurrency**: maximises throughput. Rejected in favour of bounded concurrency (`--parallel-downloads N`) to cap peak disk usage at N × ~2GB.

**Resumable Stage 2**: skip already-written traits on re-run. Rejected — Stage 2 has no downloads or heavy computation; a clean rebuild is cheaper than the deduplication logic required for safe partial writes.

## Source locator

The `file_path` column in the annotation TSV is a **source locator** — a URL or a local path. Stage 1 downloads and deletes only URLs; local paths are filtered in place and never deleted. Detection is purely syntactic (URL scheme present or absent).

## Consequences

- Stage 1 is safely resumable at trait granularity with no user intervention
- Peak disk usage is bounded to `N × source_file_size + all intermediates`
- Intermediate files are inspectable without special tooling
- Stage 2 must be rerun from scratch if interrupted, but is fast
- `discover-study` output TSV is the contract between metadata discovery and Stage 1
