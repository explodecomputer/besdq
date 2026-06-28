# Make reference completion an optional release enhancement

BESDQ ingestion first builds an Observed-Only Store Release that is immediately queryable. Reference completion is an independent optional phase that can produce a Reference-Completed Store Release by adding imputed associations against a declared reference variant set while preserving the ability to distinguish observed from imputed associations.

Published Store Releases remain immutable: enhancement produces a new release rather than mutating a published one in place. This keeps the query engine usable for both phases and allows build pipelines to stop at observed-only storage when imputation is unnecessary, too expensive, or not yet validated.

The enhanced Reference-Completed release keeps the same Store identity as the Observed-Only source release and receives a new release identity.

Observed-Only releases have implicit observed status. Reference-Completed releases must expose Association Status so queries can distinguish missing, observed, and imputed associations.

Reference completion always produces the full canonical statistic pair, Z and SE, for imputed associations. There is no Z-only completion mode in the Store contract.

The LD Reference Panel used for completion defines the Reference Variant Set for the completed release. This makes cross-store alignment depend on explicit panel identity and version rather than each Store inventing its own completion grid.

The Reference Completion Method, including the algorithm and parameters used for imputed Z and SE, is recorded once as release-level provenance. It is not repeated per association because imputed associations in a completed release share the same method contract.

Reference Completion Quality is recorded at LD-block-by-Analysis granularity rather than per association. This matches the LD-block imputation process and avoids storing noisy per-cell quality metadata.

Observed associations outside the Reference Variant Set are retained in ragged overflow rather than discarded. The dense completed grid covers the LD panel; overflow preserves source-faithful associations that cannot be placed on that grid.

For Dense Reference-Completed releases, the dense matrix axis contains only Reference Variant Set variants. Observed off-panel associations are stored exclusively in Ragged Overflow so dense axes remain identical across Stores completed with the same LD Reference Panel.

Observed-Only Dense releases remain source-faithful and do not require an LD Reference Panel or panel-defined dense axis. A later Reference-Completed release may use a different dense axis defined by the LD Reference Panel.

For Ragged Cis-and-Signals molecular Stores, reference completion is bounded to the retained cis regions and significant trans regions. These regions are completed to reference-panel variants within their existing boundaries only; singleton suggestive associations remain as observed associations and are not region-expanded.

Ragged Reference-Completed stores use the same association sequence for observed, imputed, and missing reference variants within each completed region. Each completed region contains the full slice of Reference Variant Set variants within that region boundary, plus observed off-panel variants inside the same boundary. Missing reference variants are retained as NaN statistic rows rather than omitted, keeping query semantics consistent with Dense Reference-Completed stores. Ragged stores do not create a separate imputed-ragged component.

Observed off-panel variants inside a Ragged Reference-Completed region have ordinary observed Association Status and are not labelled as a separate Query Component.
