# Include imputed associations by default in reference-completed queries

Queries against Reference-Completed Store Releases include imputed associations by default because the purpose of reference completion is to avoid missingness and replace slow query-time LD proxy lookup. Returned results must expose Association Status so callers can distinguish observed, imputed, and missing values, and query APIs should provide an observed-only mode for analyses that require source-observed associations only.

For exact and range queries against Dense Reference-Completed releases with Ragged Overflow, overflow results are included by default so observed off-panel associations do not silently disappear. Returned rows must expose which Query Component produced the result.

Top-hit queries give Dense, Ragged, and Ragged Overflow associations the same priority. Results are merged and ranked by significance rather than preferring one storage component over another; returned rows still expose Association Status and Query Component.

Both Dense and Ragged components should provide top-hit indexes using the same significance thresholds and result contract. The physical index encoding may differ by layout, but the query API exposes one merged result stream.

Imputed associations are eligible for top-hit results by default in Reference-Completed releases. The observed-only query mode applies to top-hit queries as well as exact and range queries.

Detailed Reference Completion Quality is not joined onto every returned imputed association by default. Query APIs may expose it when requested, while default result rows keep Association Status and Query Component as the minimal provenance fields.
