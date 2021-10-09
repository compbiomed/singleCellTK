# Single Cell TK
  <!-- badges: start -->
[![R-CMD-check](https://github.com/compbiomed/singleCellTK/workflows/R-CMD-check/badge.svg)](https://github.com/compbiomed/singleCellTK/actions)
[![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)
<!-- badges: end -->

The Single Cell ToolKit (SCTK) is an analysis platform that provides an **R interface to several popular single-cell RNA-sequencing (scRNAseq) data preprocessing, quality control, and visualization tools**. SCTK imports raw or filtered counts from various scRNAseq technologies and upstream tools such as 10x CellRanger, BUStools, Optimus, STARSolo, and more. By integrating several publicly available tools written in R as well as Python, SCTK performs extensive quality control measures including doublet detection and batch effect correction. Additionally, SCTK summarizes results and related visualizations in a comprehensive R markdown and/or HTML report. SCTK provides a standardized single cell analysis workflow by representing the counts data and the results using the [SingleCellExperiment](https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) R object. Furthermore, SCTK enables seamless downstream analysis by exporting data and results in flat .txt and Python Anndata formats.  

A comprehensive list of available functions is listed in the Reference section. More information about the toolkit can be found at the toolkit [homepage](https://camplab.net/sctk/).

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

If you face any difficulty in installing or using the toolkit or have identified a bug in the toolkit, please feel free to open up a [GitHub Issue](https://github.com/compbiomed/singleCellTK/issues).
