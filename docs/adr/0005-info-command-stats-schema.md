# ADR 0005 — Info command: tier stats stored in `epi`, distance-based peak fallback, hard error on old schema

## Status

Accepted

## Context

Issue #3 requests an info command that reports, per trait, how many associations came from the cis tier, how many independent significant-trans peaks were found, and how many suggestive-trans singletons were retained — alongside the total SNP count in the original source file.

Three decisions required justification.

## Decisions

### 1. Tier stats stored as columns on `epi`, not a separate table

**Alternatives considered:**
- Separate `trait_stats` table joined on `probe_idx`
- JSON blob per trait inside `epi` or `besd_meta`

**Decision:** four new columns directly on `epi`: `n_source_snps INTEGER`, `n_cis INTEGER`, `n_sig_trans_peaks INTEGER`, `n_sug_trans INTEGER`.

**Rationale:** these are scalar, per-trait, build-time attributes — the same category as `trait_var`, `trait_chr`, and `trait_bp` already on `epi`. Putting them on `epi` keeps the schema flat, makes them trivially queryable without a join, and is consistent with the existing column set. A separate table would add a join for no benefit given these are simple integers.

### 2. Distance-based peak approximation when plink2 clumping was not run

When `--ld-reference` is not supplied, plink2 clumping is skipped and all sig-trans candidates are stored as-is. Reporting a raw SNP count as "peaks" would be misleading. 

**Decision:** approximate independent peak count by counting variants where no other variant on the same chromosome is within `sig_radius` bp. This uses only stored ESI data and the `sig_radius` build parameter (stored in `besd_meta`).

**Trade-off:** the approximation can over-count in regions of high LD (two variants 600 kb apart may not be independent). The output labels this as an estimate. Users who need accurate peak counts should supply `--ld-reference` at build time.

### 3. Hard error on databases missing the new `epi` columns

**Alternatives considered:**
- Show `N/A` for missing stats with a "rebuild required" note (graceful degradation)
- Silently show only what is available

**Decision:** `--info` errors out immediately if the new columns are absent, with the message: *"database schema is outdated — please rebuild with `import-gwas-ssf`"*.

**Rationale:** partial stats are more dangerous than no stats. A user seeing `N/A` for tier breakdown may draw incorrect conclusions about what data is stored. A clear error forces a rebuild, after which the stats are complete and trustworthy. The basic counts (`n_snps`, `n_traits`) are available via other means (`--info` on the rebuilt db); there is no loss that justifies the risk of misleading partial output.

## Build parameters stored in `besd_meta`

`cis_radius`, `sig_threshold`, `sug_threshold`, `sig_radius` are written to `besd_meta` at build time. They are needed by the distance-based peak approximation and by the info display to show what filtering was applied.
