# 0006. Replace deprecated scuttle and scran functions

- **Status:** Proposed
- **Date:** 2026-10-07
- **Deciders:** TODO (maintainer approval required)

## Context

scuttle 1.22 and scran 1.40 (Bioc 3.23) deprecated most of their compute
functions in favour of scrapper and bluster. singleCellTK called many of
them directly or through scater:

- `aggregateAcrossCells` and `aggregateAcrossFeatures`
- `logNormCounts`
- `addPerCellQC`
- `buildSNNGraph`
- `pairwiseWilcox`
- `modelGeneVar`
- `quickCluster`

The deprecation warnings in examples fail R CMD check on GitHub CI, and
these functions are expected to become defunct in a later Bioconductor
release. The maintainer wants the migration in the Bioc 3.24 release.

Some replacements give identical results. Others (`modelGeneVar`,
`quickCluster`) use different algorithms, so their results change. The
replacements also don't always produce everything singleCellTK reports.
`scrapper::computeRnaQcMetrics()` lacks subset detected counts, top-N
percentages, and altExp metrics, and `pairwiseWilcox()` has no replacement
with p-values.

## Decision

We will replace each deprecated call so that singleCellTK's outputs keep
their names, shapes, and values. Where the call changes results, we change
it in a separate PR that the maintainer reviews scientifically.

- **Imports.** Add `scrapper (>= 1.6.0)`, `bluster`, and `BiocSingular`
  (for the centered SVD scran used before building SNN graphs). Remove
  `scuttle`, which has no direct uses left (scater still depends on it).
  `scran` stays, for `computeSumFactors()` (TSCAN size factors) and
  `denoisePCANumber()` (number of PCs for SoupX's quick clustering). Neither
  is deprecated.
- **Results unchanged (this PR):**
  - aggregation: scrapper, or base R `rowsum()` for TSCAN centroids;
  - `scaterlogNormCounts()`: `scrapper::normalizeRnaCounts.se()`;
  - per-cell QC: `scrapper::computeRnaQcMetrics()` for totals and detected
    counts, with internal code for the columns it lacks;
  - SNN graphs: `bluster::makeSNNGraph()`;
  - Wilcoxon DE: an internal vectorized rank-sum test.
- **Results change (separate PR):** `runModelGeneVar()` moves to scrapper's
  variance modelling, and `runSoupX()`'s clustering to a scrapper/bluster
  pipeline.
- The internal replacements live in `R/scuttleScranReplacements.R` and
  `R/runPerCellQC.R`. Each is tested against base R.

## Consequences

- No deprecation warnings from singleCellTK's own calls into scuttle or
  scran for the replaced functions.
- batchelor 1.28's `mnnCorrect()` still calls deprecated scuttle functions.
  `runMNNCorrect()` muffles only `deprecatedWarning`s raised during that
  call, with a comment saying to remove it once batchelor updates, so that
  the examples pass R CMD check. decontX 1.10 (Bioc 3.23 only) still warns.
- singleCellTK now maintains small numerical routines (QC metrics,
  Wilcoxon test) that scuttle and scran used to provide. They are short
  and tested, but they are ours to fix.
- Results of the replaced functions are identical, or within
  floating-point error (2e-15), on `scExample` and `sceBatches`.
- `runModelGeneVar()` HVG rankings and `runSoupX()`'s automatic clusters
  change. The before/after comparison is in
  `dev/plans/2026-10-07-scrapper-migration-results.md`.

## Alternatives considered

- **Keep the deprecated calls and suppress the warnings.** Rejected: they
  will become defunct, and suppressing deprecation warnings hides the next
  break.
- **Use `scrapper::scoreMarkers()` AUC for Wilcoxon DE.** Rejected: it
  has no p-values or FDR, which `runWilcox()` reports and filters on.
- **Use `scrapper::quickRnaQc.se()` for QC.** Rejected: it doesn't produce
  several columns used by QC plots, the Shiny app, and reports.
- **Use `presto` for Wilcoxon.** Rejected: it is not on CRAN or
  Bioconductor, so it can't be an Import.
