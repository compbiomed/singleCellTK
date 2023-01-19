# The Single Cell Tool Kit

[![BioC-check](https://github.com/compbiomed/singleCellTK/actions/workflows/BioC-check.yaml/badge.svg?branch=master)](https://github.com/compbiomed/singleCellTK/actions/workflows/BioC-check.yaml) [![R-CMD-check](https://github.com/compbiomed/singleCellTK/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/compbiomed/singleCellTK/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/devel/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)

The Single Cell Toolkit (SCTK) in the *singleCellTK* R package is an analysis platform that provides an **R interface to several popular single-cell RNA-sequencing (scRNAseq) data preprocessing, quality control, analysis, and visualization tools**. SCTK imports raw or filtered counts from various scRNAseq preprocessing tools such as 10x CellRanger, BUStools, Optimus, STARSolo, and more. By integrating several publicly available tools written in R or Python, SCTK can be used to perform extensive quality control including doublet detection, ambient RNA removal. SCTK integrates analysis workflows from popular tools such as Seurat and Bioconductor/OSCA into a single unified framework. Results from various workflows can be summarized and easily shared using comprehensive HTML reports. Lastly, data can be exported to Seurat or AnnData object to allow for seamless integration with other downstream analysis workflows.

![](https://camplab.net/sctk/img/interior-2.png)

## Features

SCTK offers multiple ways to analyze your scRNAseq data through the R console, the command line, and a graphical user interface (GUI) with the ability to use a large number of algorithms from both R & Python integrated within the toolkit.

-   **Interactive Analysis:** The Shiny APP allows users without programming experience to easily analyze their scRNAseq data with a GUI. Try it out at <https://sctk.bu.edu/>.

-   **Console Analysis:** Traditional analysis of scRNAseq data can be performed in the R console using wrapper functions for a multitude of tools and algorithms.

-   **Reports:** Comprehensive HTML reports developed with RMarkdown allows users to document, explore, and share their analyses.

-   **Interoperability:** Tools from both R and Python package can be seamlessly integrated within the same analysis workflow without the need for manual conversion between different data objects and file formats.

-   **Number of tools:** SCTK provides access to the largest number of tools within the same platform streamlining end-to-end analysis workflows. Curated workflows include those from Seurat, Scanpy, Scater/Scran (Bioconductor), and Celda.

## Tutorials

-   [Import and QC:](https://camplab.net/sctk/current/articles/import_data.html) The Import and QC workflows allow users to import data from multiple formats and perform comprehensive quality control and filtering.
-   ["*A la carte*" workflow:](https://camplab.net/sctk/current/articles/02_a_la_carte_workflow.html) The "A la carte" workflow lets users choose from a variety of options during each step of the analysis workflow including normalization, batch correction (optional), dimensionality reduction, 2-D embedding, and clustering.
-   [Seurat curated workflow:](https://camplab.net/sctk/current/articles/seurat_curated_workflow.html) This curated workflows recapitulates the steps for clustering and integration from the Seurat package (R).
-   [Scanpy curated workflow:](https://camplab.net/sctk/current/articles/scanpy_curated_workflow.html) This curated workflows recapitulates the steps for clustering from the Scanpy package (Python).
-   [Celda curated workflow:](https://camplab.net/sctk/current/articles/celda_curated_workflow.html) The curated Celda workflow performs matrix factorization by clustering genes into co-expression modules, cells into subpopulations, and estimating the amount of each module in each cell population.

## Installation

R package `singleCellTK` is available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html) and can be installed with the following commands:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("singleCellTK")
```

Additional information on how to install from GitHub, install Python dependencies, and for troubleshooting is available on the [Installation](https://camplab.net/sctk/current/articles/installation.html) page.

## Citation

If you use SCTK for quality control, please cite our *Nature Communication* paper

> Rui Hong, Yusuke Koga, Shruthi Bandyadka, Anastasia Leshchyk, Yichen Wang, Vidya Akavoor, Xinyun Cao, Irzam Sarfraz, Zhe Wang, Salam Alabdullatif, Frederick Jansen, Masanao Yajima, W. Evan Johnson & Joshua D. Campbell, "Comprehensive generation, visualization, and reporting of quality control metrics for single-cell RNA sequencing data," *Nature Communications*, vol. 13, no. 1688, 2022, doi: 10.1038/s41467-022-29212-9.

If you use SCTK for analysis in the Rconsole or the interactive graphical user interface, please cite our bioRxiv paper:

> Yichen Wang, Irzam Sarfraz, Rui Hong, Yusuke Koga, Vidya Akavoor, Xinyun Cao, Salam Alabdullatif, Nida Pervaiz, Syed Ali Zaib, Zhe Wang, Frederick Jansen, Masanao Yajima, W Evan Johnson, Joshua D Campbell, "Interactive analysis of single-cell data using flexible workflows with SCTK2.0", *bioRxiv*, 2022.07.13.499900; doi: <https://doi.org/10.1101/2022.07.13.499900>.

## Report Issues

If you face any difficulty in installing or have identified a bug in the toolkit, please feel free to open up an [Issue](https://github.com/compbiomed/singleCellTK/issues) on GitHub. Questions about how to best analyze your scRNA-seq data can be asked in the [Discussions](https://github.com/compbiomed/singleCellTK/discussions) page on GitHub.
