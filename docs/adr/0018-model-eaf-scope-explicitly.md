# Model EAF and INFO scope explicitly

BESDQ represents effect allele frequency with an explicit scope rather than assuming it is part of the canonical variant table. Dense Stores may store one EAF value per variant when the value is genuinely shared, or one value per association when it varies across Analyses; Ragged Stores store EAF per association because retained associations are already Analysis-specific sequences.

Imputation INFO uses the same scope model independently from EAF. EAF and INFO are optional metadata and are not used to reconstruct beta, SE, Z, or p-value in the new Store contract. Variant-scoped EAF or INFO is allowed only when the builder can establish that one value is genuinely shared; differing values are represented at association scope or omitted, not averaged.

For imputed associations, stored EAF comes from the LD Reference Panel rather than being inferred from the source study. In v1, EAF provenance is inferred from Association Status: observed associations use source EAF and imputed associations use reference-panel EAF. A separate EAF source field is deferred until mixed provenance within observed associations is required.
