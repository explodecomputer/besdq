# Use a controlled stored effect scale and loose original effect scale

BESDQ stores beta on a standardised **Stored Effect Scale** with a strict controlled vocabulary, rather than preserving every source unit as the primary storage scale. The initial vocabulary is SD Units, Log Odds, and Log Hazard; there is no `other`, `unknown`, or original-units stored scale. Unsupported scales should fail ingestion until the vocabulary is deliberately extended.

Sample size semantics remain separate from effect scale. BESDQ records whether sample size is participant, case-control, effective, or unknown, but does not require a separate effective-N method ontology at this stage.
