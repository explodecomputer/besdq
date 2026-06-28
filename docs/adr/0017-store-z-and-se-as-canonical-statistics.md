# Store Z and SE as canonical association statistics

New Dense and Ragged BESDQ stores use Z and SE as the canonical logical statistics for each retained association, with beta derived as `z × se` when needed. SE is always non-negative; effect direction is carried by signed Z. This preserves the source uncertainty directly, avoids the undefined `beta / z` case when Z is zero, and keeps reconstruction independent of allele frequency and sample size.

This supersedes ADR-0016. Beta remains queryable, but it is derived on the Store's **Stored Effect Scale** rather than stored as the primary association statistic.
