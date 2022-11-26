# Single Cell TK

[![BioC-check](https://github.com/compbiomed/singleCellTK/actions/workflows/BioC-check.yaml/badge.svg?branch=master)](https://github.com/compbiomed/singleCellTK/actions/workflows/BioC-check.yaml)
[![R-CMD-check](https://github.com/compbiomed/singleCellTK/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/compbiomed/singleCellTK/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/compbiomed/singleCellTK/branch/devel/graph/badge.svg)](https://codecov.io/gh/compbiomed/singleCellTK)

The Single Cell ToolKit (SCTK) is an analysis platform that provides an **R interface to several popular single-cell RNA-sequencing (scRNAseq) data preprocessing, quality control, analysis, and visualization tools**. SCTK imports raw or filtered counts from various scRNAseq preprocessing tools such as 10x CellRanger, BUStools, Optimus, STARSolo, and more. By integrating several publicly available tools written in R or Python, SCTK can be used to perform extensive quality control including doublet detection, ambient RNA removal. SCTK integrates analysis workflows from popular tools such as Seurat and Bioconductor/OSCA into a single unified framework. Results from various workflows can be summarized and easily shared using comprehensive HTML reports. Lastly, data can be exported to Seurat or AnnData object to allow for seamless integration with other downstream analysis workflows. More information about the toolkit can be found at the toolkit [Homepage](https://camplab.net/sctk/).

## Features

SCTK offers multiple ways to analyze your scRNAseq data both through the R console, commandline (QC) and graphical user interface (GUI) with the ability to use a large number of algorithms from both R & Python integrated within the toolkit. 

#### Interactive Analysis

The Shiny APP allows users without programming experience to easily analyze their scRNAseq data with a GUI. Try the instance at https://sctk.bu.edu/

#### Console Analysis

Traditional analysis of scRNAseq data can be performed in the R console using wrapper functions for a multitude of tools and algorithms.

#### Reports

Comprehensive HTML reports developed with RMarkdown allows users to document, explore, and share their analyses.

#### Interoperability

Tools from both R and Python can be seamlessly integrated within the same analysis workflow.

## Installation

R package `singleCellTK` is available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html). 

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("singleCellTK")
```

Detailed instruction on how to install SCTK and additional dependencies are available at our [Homepage](https://camplab.net/sctk/).

## Citation

If you use SCTK, please cite our *Nature Communication* paper

> Rui Hong, Yusuke Koga, Shruthi Bandyadka, Anastasia Leshchyk, Yichen Wang, Vidya Akavoor, Xinyun Cao, Irzam Sarfraz, Zhe Wang, Salam Alabdullatif, Frederick Jansen, Masanao Yajima, W. Evan Johnson & Joshua D. Campbell, “Comprehensive generation, visualization, and reporting of quality control metrics for single-cell RNA sequencing data,” *Nature Communications*, vol. 13, no. 1688, 2022, doi: 10.1038/s41467-022-29212-9.

## Report Issues

If you face any difficulty in installing or have identified a bug in the toolkit, please feel free to open up an [Issue](https://github.com/compbiomed/singleCellTK/issues) on GitHub. Questions about how to best analyze your scRNA-seq data can be asked in the [Discussions](https://github.com/compbiomed/singleCellTK/discussions) page on GitHub. 
