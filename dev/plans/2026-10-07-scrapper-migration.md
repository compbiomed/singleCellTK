# Migrate off deprecated scuttle/scran functions (Bioc 3.24)

## Context

scuttle 1.22 and scran 1.40 (Bioc 3.23) deprecated most of their compute
functions in favour of scrapper and bluster. The deprecation warnings in
examples make GitHub CI fail (`checking examples ... WARNING`), and the
functions will become defunct in a later release. The maintainer wants all
of it in the Bioc 3.24 release (singleCellTK 2.24). The work is split so
that changes which alter results get their own review.

The replacements were checked against the installed scrapper 1.6.3 and
bluster 1.22.0 by running old and new code side by side (research notes are
in the session record).

## PR B1: `fix/scrapper-migration` (results unchanged)

Stacked on PR A (#796), and becomes version 2.23.2. Each step is one
commit, test-first: write or extend a test that pins the current output,
watch it pass on the old code, switch the implementation, and watch it
still pass with no deprecation warning.

1. **ADR 0006.**
   - Add `scrapper`, `bluster`, and `BiocSingular` to Imports.
   - Remove the `scuttle` importFrom.
   - Drop `scuttle` and `scran` from Imports once nothing uses them. scran
     remains only if `computeSumFactors` stays, which B2 decides.
2. **Aggregation.**
   - `runClusterSummaryMetrics()`: use `scrapper::aggregateAcrossCells()`
     sums, detected counts, and cell counts for the mean and proportion
     detected. Keep the current matrix dimnames.
   - `runTSCAN` (4 sites): compute reducedDim centroids with `rowsum()`.
     Only centroids were ever used.
   - `plotSCEHeatmap()`: matrix aggregation for cells and features, rebuild
     the colData/rowData from `combinations`, and drop NA ids first (scuttle
     dropped them; scrapper keeps them).
3. **`scaterlogNormCounts()`.** Call `scrapper::normalizeRnaCounts.se()`:
   - Pass existing `sizeFactors()` when present, as scuttle reused them.
   - Use `delayed = FALSE`, so the assay stays a `dgCMatrix`.
   - Stop with an informative error for zero-count cells, as scuttle did;
     scrapper would return NaN.
   - Store the size factors in the same `sizeFactor` column.
4. **Per-cell QC** (`runPerCellQC()`, `sampleSummaryStats()`).
   - Add an internal `.perCellQCMetrics()` built on
     `scrapper::computeRnaQcMetrics()`. Fill in the columns it lacks by hand:
     `percent.top_N`, `subsets_*_sum/detected/percent`,
     `altexps_*_sum/detected/percent`, and `total`.
   - Apply `detectionLimit` by hand.
   - Column names and values stay identical.
   - Update the version recorded at `runPerCellQC.R:319`.
5. **SNN graphs** (`runScranSNN()`).
   - Use `bluster::makeSNNGraph()` for the reducedDim branches.
   - For the assay and altExp-assay branches, run
     `BiocSingular::runPCA()` first, as scran did internally, under the
     existing seed.
6. **Wilcoxon DE** (`runDEAnalysis.R:691`).
   - Add an internal vectorised two-sided Wilcoxon rank-sum test: normal
     approximation with tie and continuity correction (matching scran and
     `wilcox.test(exact = FALSE, correct = TRUE)`), `p = 1` for constant
     genes, and BH FDR. Process rows in chunks.
   - The test checks the p-values against `stats::wilcox.test()` per gene
     on a small matrix.
7. **Docs and release.**
   - Run `make docs` and keep only the relevant .Rd changes.
   - Update roxygen and vignette text that names `scran::buildSNNGraph`.
   - Add a NEWS entry ("no change to results") and bump to 2.23.2.

Verification: `make test`, `make coverage` (no drop), `make check-full`
(no deprecation WARNING from these call sites), and `make bioccheck`. The
remaining warnings, from batchelor's `mnnCorrect` and decontX 1.10 on
3.23, are upstream issues.

## PR B2: `fix/hvg-soupx-scrapper` (results change; maintainer review)

Stacked on B1, and becomes 2.23.3.

1. **`runModelGeneVar()`.**
   - Use `scrapper::modelGeneVariances()`. Map `means` → `_mean`,
     `variances` → `_totalVariance`, and `residuals` → `_bio`.
   - The rowData names and metadata stay the same.
   - The HVG ranking changes, because scrapper uses a different trend fit.
2. **`runSoupX()` clustering.**
   - Replace `scran::quickCluster()` with normalize → `chooseRnaHvgs.se()`
     (top 500) → `runPca()` → `bluster::makeSNNGraph()` →
     `igraph::cluster_walktrap()`.
   - Keep the `SoupX_cluster` name.
   - The clusters, and therefore SoupX's rho, will change.
3. **`runTSCAN()` size factors.** Switch `scran::computeSumFactors()` to
   `scuttle::computePooledFactors()` (same code), or keep scran. Decide by
   whether scran can then be dropped.
4. **Before/after report** in `dev/plans/`, on `scExample`, `sceBatches`,
   and the vignette data:
   - HVG overlap (top 500/2000) and rank correlation of bio.
   - SoupX cluster ARI and rho before/after.
   - Downstream cluster ARI after runTSNE's automatic HVG step.
5. **Vignettes.** Render the affected pkgdown articles with
   `make article` for the maintainer to review.
6. NEWS entry that names the result changes, then bump.

## PR C: `fix/shiny-celda-feature-selection` (from devel)

The app's Celda tab feature selection is broken in both branches
(`inst/shiny/server.R` ~4212–4232):

- It calls `scranModelGeneVar()`, which doesn't exist. It should be
  `runModelGeneVar(useAssay = ...)`.
- It passes `n =` to `getTopHVG()`, whose argument is `hvgNumber`. It
  also needs `useFeatureSubset = NULL`.
- It calls the deprecated `scater::logNormCounts` directly. It should use
  `scaterlogNormCounts()`.

The maintainer explicitly asked for this app change. It is verified with a
screenshot of the running app (`make app`), NEWS, and a z bump.

## Version order

PR A 2.23.1 → B1 2.23.2 → B2 2.23.3 → C 2.23.4. If C merges before B2, the
numbers are renumbered at merge time. All must be pushed to Bioc before
the 3.24 freeze.
