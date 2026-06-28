# BESDQ Domain Glossary

## GCST Accession

The EBI GWAS Catalog identifier for a single GWAS-SSF file. In BESDQ's domain, each GCST accession identifies exactly one **Analysis**. Multiple GCST accessions may be grouped into a **Study Batch**.

## Analysis

A particular statistical analysis of one **Trait**, producing associations between that trait and a set of variants. The same trait may have multiple analyses that differ by cohort, ancestry, model, sample subset, or meta-analysis.

## Study Batch

A set of GCST accessions sharing a PubMed ID. The natural unit of bulk ingestion. Intermediate files from Stage 1 are grouped in a subdirectory named by PMID (or a user-supplied project name).

## BESDQ Store

A self-contained logical distribution unit containing one or more **Analyses** and everything required to interpret and query them. A Store has a stable identity, follows the common BESDQ contract, and has a **Primary Storage Layout**.

## Store Release

An immutable, self-identifying published version of a **BESDQ Store**. Each release records its Store identity, release identity, format version, creation time, and provenance so that downloaded or mirrored copies remain interpretable without a catalogue service.

## Format Version

The version of the BESDQ storage contract required to interpret a **Store Release**. It describes compatibility of the representation, not the version of the data contained in the Store.

## Primary Storage Layout

The main physical organisation of associations within a **BESDQ Store**. Dense and Ragged are alternative primary layouts behind the same metadata, identity, validation, and query concepts; layout is independent of **Association Coverage**.

## Dense Layout

A **Primary Storage Layout** in which Analyses share a variant axis and associations occupy cells in a variant-by-Analysis matrix. An explicit missing value represents an unavailable association.

For Reference-Completed Dense stores, the dense variant axis is the **Reference Variant Set**; observed off-panel variants are retained in **Ragged Overflow**.

## Ragged Layout

A **Primary Storage Layout** in which each Analysis has its own sequence of retained associations referencing a Store-wide variant table. Different Analyses may retain different variant sets.

In Reference-Completed Ragged stores, observed, imputed, and missing reference variants for a completed region belong to the same Analysis association sequence and are distinguished by **Association Status**.

## Ragged Overflow

An optional ragged component of a Dense Reference-Completed **Store Release** that preserves observed associations outside the **Reference Variant Set**. The dense component remains the primary query grid.

## Query Component

The part of a multi-component Store Release from which a returned association was read, such as Dense Grid or Ragged Overflow. Off-panel observed variants inside a Ragged primary layout's completed region are ordinary observed associations, not a separate Query Component.

## Association Coverage

The guarantee a **Store Release** makes about which source associations it retains, independently of its **Storage Layout**. Full Coverage retains every usable source association after normalisation and quality control; Cis-and-Signals Coverage retains complete cis regions plus selected significant and suggestive trans associations.

## Completion State

Whether a **Store Release** contains only source-observed associations or has also been completed against a reference variant set. Observed-Only and Reference-Completed releases share the same query concepts but differ in whether imputed associations may be present.

## LD Reference Panel

The declared ancestry-specific LD resource used for reference completion. It defines the **Reference Variant Set** and provides the LD information used to infer imputed associations.

## Reference Variant Set

The canonical variant set defined by an **LD Reference Panel** for a **Reference-Completed** release. Reference completion attempts to provide associations on this set, subject to missingness where imputation fails or is out of scope.

## Reference Panel Membership

Store variant metadata indicating whether a variant belongs to the **Reference Variant Set** for a Reference-Completed release. It distinguishes reference-panel variants from observed off-panel variants.

## Reference Completion Region

A genomic interval within which a Ragged Cis-and-Signals Store attempts reference completion. A completed region contains every **Reference Variant Set** variant inside its boundary plus observed off-panel variants in the same boundary; singleton suggestive associations are not expanded.

## Reference Completion Method

The release-level algorithm and parameters used to infer imputed Z and SE values from observed associations and the **LD Reference Panel**. It is recorded as provenance for a Reference-Completed release rather than repeated per association.

## Reference Completion Quality

A quality summary for imputed associations at LD-block-by-Analysis granularity. It describes confidence in a block of imputed values rather than attaching separate quality metadata to every association.

## Observed Association

An association whose statistics come from the source dataset after BESDQ normalisation. Observed associations are the authoritative basis for a Store Release.

## Imputed Association

An association whose Z and SE were inferred during reference completion rather than reported by the source dataset. Imputed associations must remain distinguishable from observed associations.

## Association Status

The origin state of an association in a Reference-Completed **Store Release**: Missing, Observed, or Imputed. Missing means the association was not observed and reference completion did not impute it.

## Observed-Only Query

A query mode that excludes **Imputed Associations** and returns only source-observed results. Reference-Completed releases include imputed associations by default, but callers may request observed-only results.

## Top-Hit Query

A query that returns associations ranked by statistical significance. Dense, Ragged, Ragged Overflow, observed, and imputed results have equal priority; ranking is determined by significance, not by storage component or association status.

## Top-Hit Index

A layout-specific acceleration structure that supports **Top-Hit Queries** using the Store's shared significance thresholds and result contract. Dense and Ragged components may encode the index differently but expose the same query semantics.

## Two-Stage Build

The pipeline for constructing a BESDQ index from EBI GWAS-SSF files at scale:

- **Stage 1** (`build-stage1`): for each Trait in a discovery TSV, download the GWAS-SSF file, extract cis/trans regions into a TSV+YAML intermediate file pair, delete the original download. Resumable: a trait is considered complete when both `GCST*.tsv.gz` and `GCST*.yaml` exist (written atomically). Failed traits produce `GCST*.failed.yaml` and are retried on re-run. Supports bounded-concurrency download via `--parallel-downloads N`.
- **Stage 2** (`build-stage2`): consume all completed intermediates in a PMID subdirectory and collate into the final SQLite index. Always a clean rebuild (no resumability needed — no downloads, fast).

## Intermediate File Pair

The on-disk representation of a single Trait's filtered associations produced by Stage 1. Two files per trait:
- `GCST<id>.tsv.gz` — filtered rows (cis + trans), human-readable gzipped TSV
- `GCST<id>.yaml` — per-trait metadata (trait_chr, trait_bp, trait_var, tier counts, etc.) needed by Stage 2

Written atomically via `.tmp` rename. Stored under `<workdir>/<pmid>/`.

## Source Locator

The value in the `file_path` column of the annotation TSV passed to `build-stage1`. Resolved at runtime as:

- **URL** (starts with `http://`, `https://`, or `ftp://`): file is downloaded to a temp location, filtered, then **deleted** after the intermediate is written.
- **Local path**: file is filtered in place and **never deleted** — the user owns it.

The distinction is purely syntactic (URL scheme present or absent). No other behaviour differs between the two cases.

## Discovery Script

The `discover-study` command. Given a PubMed ID, queries the EBI GWAS Catalog API to enumerate all GCST accessions and their FTP URLs. Optionally parses gene symbols from trait description fields (`--parse-gene` flag, default off) and maps them to chr/bp via a local Ensembl gene annotation file. Emits an annotation TSV consumed by `build-stage1`.

## Trait

A measured or derived outcome, such as disease status, LDL cholesterol, gene expression, methylation, or protein abundance. **Phenotype** is a synonym, but BESDQ uses Trait consistently; one Trait may be the subject of multiple **Analyses**.

**EPI columns by class:**

| Class | Columns | Behaviour when absent |
|---|---|---|
| Mandatory | `trait_id TEXT`, `trait_name TEXT` | Import fails |
| Optional functional | `trait_var REAL`, `trait_chr TEXT`, `trait_bp INTEGER` | Affects reconstruction / cis-window; warned at build time |
| Optional non-consequential | `gene TEXT`, `context TEXT` | NULL; no effect on queries or reconstruction |

`gene`: human-readable gene symbol (e.g. IL10). `context`: free-text descriptor covering tissue, treatment, time-point, sex, or any other experimental condition that distinguishes traits within a study (e.g. "PBMC_Bbmix_baseline").

`orientation` is dropped — it was carried from the BESD `.epi` convention but is unused in any query or output.

Study-level metadata (publication, year, ancestry, max sample size, study type, source URL) is stored as a single JSON blob in the `besd_meta` table under the key `study_metadata`. This is dataset-wide, not per-trait.

> **Probe** in the BESD binary format is the historical synonym for Trait. The `.epi` file, `probe_data` table, and binary block structure still use "probe" internally.

## Variant Identity and Allele Convention

The canonical within-Store variant key is the ALID `chr:pos:A1:A2`, where **A1 is always the alphabetically first allele** and A2 is the other. Cross-Store identity is the pair (**Reference Assembly**, ALID); the rsid is an alias rather than primary identity, and variants without an rsid remain valid.

Every variant-Analysis association is normalised to this orientation: if the source effect allele is not alphabetically first, the alleles are swapped and signed statistics are negated.

Multiple canonical variants may share a genomic position when their allele pairs differ. A multiallelic source record is decomposed when it provides unambiguous per-allele association statistics; an ambiguous record does not define a canonical association.

Alleles are trimmed and left-aligned against the Store's **Reference Assembly** before identity is assigned. A long-allele ALID may use a deterministic hash as its compact representation, but the complete normalised alleles remain authoritative and are retained by the Store.

## Store Variant Table

The Store-wide union of canonical variants referenced by its **Analyses**. Each variant occurs once in a **Store Release**, independently of how many Analyses report it.

## Variant Index

A compact Store-local reference to a row in the **Store Variant Table**. It has no identity or stability guarantee outside its Store Release; cross-Store matching uses canonical variant identity instead.

## Reference Assembly

The genome assembly to which every variant coordinate in a **Store Release** refers, such as GRCh37 or GRCh38. Each Store Release declares exactly one Reference Assembly; coordinate conversion occurs before ingestion and is part of provenance.

## Statistics Encoding

The logical association statistics stored by a BESDQ Store. New Dense and Ragged Stores store Z and SE for each retained association; beta is derived from those statistics on the same **Stored Effect Scale**.

## Z-Score

`z = beta / se`. Z is signed, carries effect direction, and is invariant to simple rescaling of beta and SE.

## Standard Error

The non-negative uncertainty of beta on the same **Stored Effect Scale**. New Dense and Ragged Stores store SE directly rather than reconstructing it from allele frequency, sample size, or beta divided by Z.

## Reconstruction

The process of deriving query output from stored association statistics. New Dense and Ragged Stores return Z and SE directly and derive beta as Z multiplied by SE when beta is requested.

## Sample Size

The source-reported number or effective number of participants contributing to an **Analysis** or to one of its variant associations. Its meaning is described by **Sample Size Kind** and its granularity by **Sample Size Scope**; BESDQ does not present a value inferred from effect statistics as an observed participant count.

## Sample Size Kind

The interpretation of a **Sample Size**: Participants, Case-Control, Effective, or Unknown. Case-Control distinguishes case and control counts; Effective is explicitly not an observed participant count.

## Sample Size Scope

Where a **Sample Size** applies: Analysis when one value applies throughout, Variant when values may differ between associations, or None when sample size is unknown.

## Stored Effect Scale

The controlled scale in which an **Analysis** stores beta and SE: SD Units for continuous traits, Log Odds for binary traits, or Log Hazard for survival traits.

## Original Effect Scale

The source-reported or source-measurement scale of beta before BESDQ normalisation, such as kg/m², mmol/L, log-odds, or source-specific free text. BESDQ records this as provenance rather than forcing it into a strict ontology.

## Phenotype Standard Deviation

The Analysis-level scale factor used to convert linear continuous effects from **Original Effect Scale** to phenotype standard deviation units. Its value and provenance are part of Analysis metadata when the conversion is meaningful.

## Trait Variance (trait_var)

Variance of the phenotype (e.g. gene expression level) for a trait. Per-trait quantity stored in `epi.trait_var REAL`. Used to compute `sd_y = sqrt(trait_var)` for SD-unit output.

**Sources, in precedence order:**
1. User-supplied in annotation TSV `trait_var` column (exact)
2. Auto-estimated at build time from cis SNPs using the median formula:
   `trait_var = median( se_i² × n × 2 × eaf_i × (1 − eaf_i) )` over all cis SNPs
3. 1.0 if no cis SNPs and not user-supplied — output is returned in original units with a warning

Auto-estimation is robust to per-SNP missing data (median estimator) and approximately correct under covariate adjustment (effective-n error ≈ k/n, typically < 1%). Previously named `var_y`.

## Output Scale

All query methods return beta and SE in **SD units** (standard deviation of the trait) by default, making effect sizes directly comparable across studies and traits. Pass `--original-scale` to return beta and SE in original study units.

SD-unit conversion: `beta_sd = beta_orig / sd_y`, `se_sd = se_orig / sd_y`, where `sd_y = sqrt(trait_var)` from EPI. The z-score and p-value are invariant — they are identical in both scales.

## Effect Allele Frequency

The frequency of the effect allele associated with an association. For observed associations EAF comes from the source dataset when available; for imputed associations EAF comes from the **LD Reference Panel** when stored.

## EAF Scope

Where an **Effect Allele Frequency** value applies: Variant when one value is shared for a Store variant, Association when values may differ by Analysis-variant association, or Absent when EAF is not stored.

## Imputation INFO

The source-reported imputation quality or information score for a variant association. INFO is optional association metadata and is not required to reconstruct beta, SE, Z, or p-value in new BESDQ Stores.

## INFO Scope

Where an **Imputation INFO** value applies: Variant when one value is shared for a Store variant, Association when values may differ by Analysis-variant association, or Absent when INFO is not stored.

## Significance Mask

The filter applied during import from text files (e.g. GWAS-SSF) to determine which SNP-trait associations are stored. Three tiers:

| Tier | Condition | What is stored |
|---|---|---|
| Cis | SNP within cis-radius of `trait_chr`/`trait_bp` | All variants unconditionally |
| Significant trans | p < 5×10⁻⁸ outside cis | All variants within sig-radius of each independent lead SNP |
| Suggestive trans | 5×10⁻⁸ ≤ p < 1×10⁻⁴ outside cis | That variant only |
| Below suggestive | p ≥ 1×10⁻⁴ outside cis | Dropped |

**Cis-radius**: the distance from `trait_bp` within which all variants are retained unconditionally. Default 1,000,000 bp (±1 Mb). Requires `trait_chr` and `trait_bp`; if absent, cis tier is skipped and only trans tiers apply.

**Sig-radius**: the distance from a significant trans lead SNP within which all variants are retained. Default 500,000 bp (±500 kb).

Independent significant trans peaks are identified by LD clumping with plink2 (must be on PATH) and a user-supplied plink2-format LD reference panel (`--pfile` prefix). Clumping defaults: r²=0.01, window=10,000 kb. Thresholds and radii are configurable at import time. When plink2 clumping is not run, peak count is approximated by counting variants more than sig-radius apart as independent.

**Peak count**: the number of LD-independent signals within a tier. Distinct from SNP count (which counts all stored variants including neighbourhood expansion around peaks).
- Cis tier: no peak concept — all variants stored unconditionally; reported as raw SNP count.
- Significant trans tier: number of independent peaks (lead SNPs from clumping, or distance-based approximation). The total SNPs stored equals all variants within sig-radius of each peak.
- Suggestive trans tier: each stored variant is already a singleton (no neighbourhood expansion), so peak count equals SNP count.

BESD files imported via the legacy path arrive pre-filtered and are stored verbatim — the Lean Index applies no additional filters to them.

## Dense Batch Format — Design Decisions

These decisions apply to the dense Zarr-based batch format (`.besdz` stores) and are format-agnostic in the sense that they would carry over to any successor format.

### Variant Identity

Same ALID convention as the Lean Index: canonical key is `chr:pos:A1:A2` where **A1 is alphabetically first**. The sign of beta and z-score must match — if the source effect allele is not alphabetically first, alleles are swapped and the sign is negated at build time. This is identical to the besdq sparse convention and to pleioDB.

### Variant Table

Built as the **union** of all variants across all traits in a batch. Variants absent from a given trait are encoded as float16 NaN in the stats arrays (no separate missingness mask). This accommodates traits with MAC filters or rare-disease exclusions, which can have substantial missingness relative to the batch-wide variant list.

EAF follows the Store's **EAF Scope**. Dense stores may use Variant or Association scope; Ragged stores use Association scope because each Analysis has its own retained association sequence.

### Stats Arrays

Two separate Zarr arrays per Dense store: `se` and `zscore`, both shape `(N_variants, N_analyses)`. Beta is not stored — derived at query time from Z and SE. P-value is derived from Z. Separate arrays allow queries that only need Z to skip loading SE.

### Effect Scale

Beta is derived on the **Stored Effect Scale**. For continuous traits this is normally phenotype standard deviation units; for binary and survival traits this is normally the model-native scale such as log-odds or log-hazard. The **Original Effect Scale** is retained as provenance rather than used as the primary storage scale.

### Sample Size N

Scalar per trait, extracted from VCF header at build time. Not stored per-variant.

### Prototype Scope

The prototype is implemented as **standalone scripts or a notebook**, with no integration into the `besdq` package. Its sole purpose is to validate storage size, query latency, and MR precision. Integration into the besdq CLI and package happens only if the benchmark confirms the format is viable. The pleioDB project serves a similar purpose to the dense batch format; its conventions (ALID, effect scale handling) are the relevant prior art to port, and it may become redundant once the dense format is production-ready.

### Store Suffix

`.besdz` (provisional). The `z` signals Zarr. The `besd` prefix may change if the project renames away from the BESD lineage.

### Significance Index (Dense Format)

A query-acceleration structure stored **in-place inside a `.besdz` store**, alongside `se` and `zscore`. Not to be confused with the Lean Index [Significance Mask](#significance-mask), which is an import-time filter that determines what gets stored; the Significance Index is a pre-computed lookup table over already-stored data.

Three p-value thresholds: 5×10⁻⁸ (genome-wide significant), 5×10⁻⁶, 5×10⁻⁴.

**Storage layout** (6 Zarr arrays total, one pair per threshold):
- `sig_5e8` — flat uint32 array of variant row indices, sorted per trait and concatenated across all N_traits traits
- `sig_5e8_offsets` — int64 array of shape `(N_traits + 1,)`; trait t's hits are `sig_5e8[offsets[t]:offsets[t+1]]`
- Same pattern for `sig_5e6` and `sig_5e4`

P-values are derived from stored Z-scores: `p = 2 × (1 − Φ(|Z|))`. Variants with NaN Z-score (untested) are excluded from all indexes.

**Tophits query** with index (vs naive full-column scan):
1. Load `sig_{key}_offsets[t:t+2]` — two int64 values
2. Slice `sig_{key}[start:end]` — the hit row indices for trait t
3. Fancy-index `se` and `zscore` at those rows — only the chunks containing hits are read

The index is built by a standalone script (`dense_11_build_sig_index.py`) and can be updated independently of the stats arrays.

## Build Modes (Lean Index)

| Dataset type | Storage mode | SE source | sd_y source |
|---|---|---|---|
| GWAS-SSF (any source with SE) | VectorN | Direct from file (float16) | Auto-estimated from cis SNPs; user override via `trait_var` |
| Single-cohort BESD, no AF | ScalarN | Reconstructed from n, AF (derived), trait_var | User-supplied `trait_var` or 1.0 |
| Single-cohort BESD, with AF | ScalarN | Reconstructed from n, AF, trait_var | User-supplied `trait_var` or 1.0 |
| Meta-analysis BESD, with AF | ScalarN (per-pair n derived at query time) | Reconstructed | User-supplied `trait_var` or 1.0 |
| Meta-analysis BESD, no AF | — | Unresolvable | **Error at build time** |
