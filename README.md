# Single Cell TK
  <!-- badges: start -->
[![R-CMD-check](https://github.com/compbiomed/singleCellTK/workflows/R-CMD-check/badge.svg)](https://github.com/compbiomed/singleCellTK/actions)
[![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/devel/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)
<!-- badges: end -->

The Single Cell ToolKit (SCTK) is an analysis platform that provides an **R interface to several popular single-cell RNA-sequencing (scRNAseq) data preprocessing, quality control, analysis, and visualization tools**. SCTK imports raw or filtered counts from various scRNAseq preprocessing tools such as 10x CellRanger, BUStools, Optimus, STARSolo, and more. By integrating several publicly available tools written in R or Python, SCTK can be used to perform extensive quality control including doublet detection, ambient RNA removal, and batch effect correction. SCTK integrates analysis workflows from popular tools such as Seurat and Bioconductor/OSCA into a single unified framework. Results from various workflows can be summarized and easily shared using comprehensive HTML reports. Lastly, data can be exported to Seurat or AnnData object to allow for seamless integration with other downstream analysis workflows. More information about the toolkit can be found at the toolkit [homepage](https://camplab.net/sctk/).

## Installation

Detailed instruction on how to install SCTK and additional dependencies are available at our homepage:
https://camplab.net/sctk/

## Features
SCTK offers multiple ways to analyze your scRNAseq data both through the R console, commandline (QC) and graphical user interface (GUI) with the ability to use a large number of algorithms from both R & Python integrated within the toolkit. 
#### Console Analysis
Traditional analysis of scRNAseq data can be performed in the R console using wrapper functions for a multitude of tools and algorithms.
#### Interactive Analysis
The Shiny APP allows users without programming experience to easily analyze their scRNAseq data with a GUI.
#### Reports
Comprehensive HTML reports developed with RMarkdown allows users to document, explore, and share their analyses.
#### Interoperability
Tools from both R and Python can be seamlessly integrated within the same analysis workflow.

## Report Issues

If you face any difficulty in installing or have identified a bug in the toolkit, please feel free to open up an [Issue](https://github.com/compbiomed/singleCellTK/issues) on GitHub. Questions about how to best analyze your scRNA-seq data can be asked in the [Discussions](https://github.com/compbiomed/singleCellTK/discussions) page on GitHub. 
