# Single Cell TK

[![Travis build status](https://travis-ci.org/compbiomed/singleCellTK.svg?branch=master)](https://travis-ci.org/compbiomed/singleCellTK)
[![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)
[![BioC status](https://www.bioconductor.org/shields/build/release/bioc/singleCellTK.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/singleCellTK)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

## Installation

### Release Version

You can download the release version of the Single Cell Toolkit in
[Bioconductor v3.10](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html):

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("singleCellTK")
```

### Devel Version

You can download the development version of the Single Cell Toolkit in
[Bioconductor v3.11](https://bioconductor.org/packages/devel/bioc/html/singleCellTK.html)
or from this repository:

```r
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK")
```

### R 3.4 Version

If you are still running an earlier version of R than 3.5, you can install
the following version from this repository:

```r
# install.packages("devtools")
devtools::install_github("compbiomed/singleCellTK", ref="r_3_4")
```

#### Troubleshooting Installation

For the majority of users, the commands above will install the latest version
of the singleCellTK without any errors. Rarely, you may encounter an error due
to previously installed versions of some packages that are required for the
singleCellTK. If you encounter an error during installation, use the commands
below to check the version of Bioconductor that is installed:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::version()
```

If the version number is not 3.6 or higher, you must upgrade Bioconductor to
install the toolkit:

```r
BiocManager::install()
```

After you install Bioconductor 3.6 or higher, you should be able to install the
toolkit using `devtools::install_github("compbiomed/singleCellTK")`. If you
still encounter an error, ensure your Bioconductor packages are up to date by
running the following command.

```r
BiocManager::valid()
```

If the command above does not return `TRUE`, run the following command to
update your R packages:

```r
BiocManager::install()
```

Then, try to install the toolkit again:

```r
devtools::install_github("compbiomed/singleCellTK")
```

If you still encounter an error, please [contact us](mailto:dfj@bu.edu) and
we'd be happy to help.

## QC Outputs
There are several available QC algorithms that are implemented within singleCellTK as wrapper functions, which will be stored as `colData` within the output `singleCellExperiment` object. These are the currently available QC outputs:

### General metrics

| Output name | Description | Package |
| --- | --- | --- |
| sum | Total transcript counts in cell | scater |
| detected | Total genes detected in cell | scater |
| percent_top | Numeric value, the percentage of counts assigned to the percent_topage of most highly expressed genes. Each column of the matrix corresponds to an entry of the sorted percent_top, in increasing order | scater |
| subsets_mito_sum | Number of total mitochonrial transcript counts per cell | scater |
| subsets_mito_detected | Number of mitochondrial genes detected per cell | scater |
| subsets_mito_percent | Percentage of mitochondial transcript counts out of total gene counts | scater |


### Metrics on Droplet matrix

| Output name | Description | Package |
| --- | --- | --- |
| dropletUtils_emptyDrops_total | Integer, spicifies the total UMI count for each barcode | dropletUtils |
| dropletUtils_emptyDrops_pvalue | Numeric, the Monte Carlo p-value under the null model | dropletUtils |
| dropletUtils_emptyDrops_logprob | Numeric, the barcode's count log-probability of a vector under the null model | dropletUtils |
| dropletUtils_emptyDrops_fdr | Numeric, false discovery rate. Suggested fdr cut-off is 1% | dropletUtils |
| dropletUtils_emptyDrops_limited | Logical, indicates if a lower p-value could be obtained by increasing niters, a number of iterations for Monte Carlo p-value calculations | dropletUtils |
| dropletUtils_BarcodeRank_Knee | Numeric, specifies total count at the knee point | dropletUtils |
| dropletUtils_BarcodeRank_Inflection | Numeric,  specifies total count at the inflection point | dropletUtils |

### Metrics for doublet detection

| Output name | Description | Package |
| --- | --- | --- |
| doubletFinder_doublet_score | Numeric value that determines how likely a cell in the counts matrix is a doublet using artificially generated doublets | doubletFinder |
| doubletFinder_doublet_label | Whether the cell is deemed a doublet or not by the algorithm. Will be "Singlet" or "Doublet" | doubletFinder |
| scds_cxds_score | Numeric value that determines how likely a cell is a doublet, based on co-expression of gene pairs | scds |
| scds_bcds_score | Numeric value that determines how likely a cell is a doublet, using artificially generated doublets | scds |
| scds_hybrid_score | Numeric value that determines how likely a cell is a doublet, uses both cxds and bcds algorithm | scds |
| scran_doubletCells_Score | Numeric value that determines how likely a cell in the counts matrix is a doublet | scran |
| scrublet_score | Numeric value that determines how likely a cell in the counts matrix is a doublet | scrublet |
| scrublet_call | Whether the cell is deemed a doublet or not by the algorithm. Will be  | scrublet |

### Metrics for ambient RNA contamination

| Output name | Description | Package |
| --- | --- | --- |
| decontX_Contamination | Probability of contamination determined by decontX | celda |
| decontX_Clusters | Clusters determined by Celda, a clustering algorithm that runs in the background of decontX | celda |

## Develop singleCellTK

To contribute to singleCellTK, follow these steps:

__Note__: Development of the singleCellTK is done using the latest version of R.

1. Fork the repo using the "Fork" button above.
2. Download a local copy of your forked repository "```git clone https://github.com/{username}/singleCellTK.git```"
3. Open Rstudio
4. Go to "File" -> "New Project" -> "Existing Directory" and select your git repository directory

You can then make your changes and test your code using the Rstudio build tools.
There is a lot of information about building packages available here: http://r-pkgs.had.co.nz/.

Information about building shiny packages is available here: http://shiny.rstudio.com/tutorial/.

When you are ready to upload your changes, commit them locally, push them to your
forked repo, and make a pull request to the compbiomed repository.

Report bugs and request features on our [GitHub issue tracker](https://github.com/compbiomed/singleCellTK/issues).

Join us on [slack](https://compbiomed.slack.com/)!
