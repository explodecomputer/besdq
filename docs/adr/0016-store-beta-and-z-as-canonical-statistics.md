---
status: superseded by ADR-0017
---

# Store beta and Z as canonical association statistics

New Dense and Ragged BESDQ stores use beta and Z as the canonical logical statistics for each retained association, with SE derived from `abs(beta / z)` when needed. This aligns dense and ragged layouts around the same query semantics and avoids using allele frequency or sample size to reconstruct effect uncertainty.

Older SQLite sparse indexes may still store Z plus SE or scalar-N reconstruction data, but those encodings are implementation history rather than the target BESDQ Store contract.
