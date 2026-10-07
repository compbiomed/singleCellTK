# Before/after results: modelGeneVar and SoupX clustering on scrapper

These are results for the maintainer's scientific review of PR B2
(`fix/hvg-soupx-scrapper`). They compare the code at
`fix/scrapper-migration` (scran's `modelGeneVar` and `quickCluster`)
with B2.

Environment: R 4.6.1, Bioconductor 3.23 packages, scrapper 1.6.3,
bluster 1.22.0, scran 1.40.0.

Data:

- `scExample` (PBMC, 195 cells without empty droplets).
- `sceBatches`.
- `scRNAseq::ZeiselBrainData()`: 3,005 mouse brain cells with published
  `level1class` cell types.

## runModelGeneVar()

`scran::modelGeneVar()` is replaced by `scrapper::modelGeneVariances()`.

- The means and total variances are identical.
- The biological component is now the residual from scrapper's LOWESS
  trend on quarter-root variances (`min.mean = 0.1`). Previously it came
  from scran's weighted loess fit.

| Data | Genes | Spearman (bio) | Genes with bio > 0 | Top-500 overlap | Top-2000 overlap |
|---|---|---|---|---|---|
| PBMC | 200 | 0.977 | 97 → 95 | 95.9% | 95.9% |
| sceBatches | 100 | 0.888 | 48 → 54 | 89.6% | 89.6% |
| Zeisel | 20,006 | 0.941 | 9,668 → 9,525 | 97.2% | 95.5% |

PBMC and sceBatches have fewer genes than 500 with bio > 0, so their
top-500 and top-2000 lists are all such genes.

Downstream on Zeisel, using 2000 HVGs → `scaterPCA` (20 PCs) →
`runScranSNN` (walktrap):

- Clusters go from 20 to 19, with an ARI of 0.777 between old and new
  labels.
- Agreement with the published cell types is 0.410 before and 0.421 after
  (ARI).

## runSoupX() automatic clustering

`scran::quickCluster(method = "igraph")` is replaced by
`.quickClusterRNA()`. It follows the same steps, using scrapper for
normalization, variance modelling, HVGs, and PCA, and bluster for the SNN
graph. It keeps scran's `denoisePCANumber()` to choose the number of PCs,
the walktrap clustering, and the merging of clusters smaller than 100
cells.

| Data | Clusters | ARI old vs new | ARI vs cell type (old → new) |
|---|---|---|---|
| Zeisel | 14 → 13 | 0.716 | 0.437 → 0.450 |
| PBMC (195 cells) | 1 → 1 | – | – |

SoupX on Zeisel (`runSoupX()`, no background, automatic clusters):

| | Before | After |
|---|---|---|
| Estimated rho | 0.013 | 0.010 |
| rho interval (FWHM) | 0.010–0.034 | 0.010–0.036 |
| Total corrected counts | 44,353,168 | 44,487,980 (+0.3%) |

The corrected counts correlate at 1.00000; the largest change in any one
entry is 119 counts. The new rho lies inside the old estimate's interval.

## Not compared

- `runSoupX()` fails on the 195-cell PBMC subset with both versions. The
  automatic clustering merges clusters below 100 cells, which leaves one
  cluster, and SoupX needs at least two (#797). The error handler then
  hides or re-raises the failure depending on the caller (#798).
- Articles were not re-rendered here. The maintainer reviews the pkgdown
  articles that use `runModelGeneVar()`: feature_selection,
  dimensionality_reduction, 2d_embedding, 02_a_la_carte_workflow,
  celda_curated_workflow, differential_expression, find_marker, and
  trajectoryAnalysis.
