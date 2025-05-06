Changes in Version 2.18.0 (2025-04-15)
================================================================================
* Updated call to msigdbr to work with newer version
* Updated version to match Bioconductor 3.21
* Fixed reproducibility with scDblFinder

Changes in Version 2.16.1 (2024-02-12)
================================================================================
* Fixed decontX/soupX webapp error
* Fixed error in reportCellQC 
* Fixed additional bugs with Seurat V5 integration

Changes in Version 2.16.0 (2024-10-29)
================================================================================
* Updated version to match Bioconductor 3.20

Changes in Version 2.14.0 (2024-05-03)
================================================================================
* Updated version to match Bioconductor 3.19
* Update runGSVA function to work with newer GSVA package

Changes in Version 2.12.2 (2024-01-28)
================================================================================
* Added support for Seurat V5

Changes in Version 2.12.1 (2024-01-10)
================================================================================
* Updates to documentation
* Fixes to runTSCAN and plotSeurat Genes
* Added support for flat file import into SCTK-QC
* Fixed directory issue in importCellRanger
* Added Bubble plot to Shiny GUI
* Updated Dockerfile

Changes in Version 2.12.0 (2023-10-24)
================================================================================
* Updated version to match Bioconductor 3.18

Changes in Version 2.10.1 (2023-07-26)
================================================================================
* Added function for bubble plot
* In SCTK-QC pipeline, added support for batch processing multiple inputs
* In SCTK-QC pipeline, added support for importing and exporting AnnData objects
* In SCTK-QC pipeline, fixed a bug causing YAML output files to be empty
* Update the SCTK-QC tutorial 
* Fixed bug in combineSCE causing it to create multiple copies of row or column data

Changes in Version 2.10.0 (2023-04-25)
================================================================================
* Updated version to match Bioconductor 3.17

Changes in Version 2.8.1 (2023-03-10)
================================================================================
* Added scanpy wrapper functions for use from console
* Added scanpy UI curated workflow
* Integrated scanpy to a la carte workflow
* Fixed a bug in importing fluidigm dataset
* Updated downloading features in the Shiny app
* Added error checking around Enrichr functions
* Minor tweaks to plot defaults

Changes in Version 2.8.0 (2022-12-19)
================================================================================
* Updated version to match Bioconductor version

Changes in Version 2.7.3 (2022-10-25)
================================================================================
* Fixed bugs related to dependency updates

Changes in Version 2.7.2 (2022-10-19)
================================================================================
* Deprecated `findMarkerDiffExp()`, `findMarkerTopTable()` and `plotMarkerDiffExp()`, which are replaced by `runFindMarker()`, `getFindMarkerTopTable()` and `plotFindMarkerHeatmap()`, respectively
* Added `useReducedDim`, `detectThresh` arguments for find marker functions
* Deprecated `getUMAP()` and `getTSNE()`, which are replaced by `runUMAP()` and `runTSNE()`, respectively
* Added `runQuickUMAP()` and `runQuickTSNE()` functions which directly compute the proper embedding from raw counts matrices with a simplified argument set
* Added arguments `aggregateRow` and `aggregateCol` to `plotSCEHeatmap()`
* Updated output metadata structure of QC functions, as well as `combineSCE()` which merges the new structure properly
* Refined batch correction function set
* Fixed bugs related to UI and console functions

Changes in Version 2.7.1 (2022-06-29)
================================================================================
* Refactored scaling related parts of the workflow, including variable feature
detection, dimension reduction and 2D embedding
* Redesigned UI landing page and UI running prompt
* Added marker table module across UI
* Added more unit tests
* Fixed bugs in TSCAN UI
* Other minor bug fixes

Changes in Version 2.6.0 (2022-04-28)
================================================================================
* Updated version to match Bioconductor

Changes in Version 2.5.2 (2022-04-23)
================================================================================
* Added Seurat report functions
* Added TSCAN trajectory analysis functions
* Refactored EnrichR wrapper function (runEnrichR)
* Added new cut-offs for DE functions
* Other refactors and bug fixes

Changes in Version 2.5.1 (2022-03-31)
================================================================================
* Added SoupX method for decontamination (runSoupX)
* Added useReducedDim parameter for DE analysis and Heatmap
* Added Differential Abundance section to the tutorials
* Fixed Mitochondrial gene list
* Other refactors and bug fixes

Changes in Version 2.4.1 (2021-12-22)
================================================================================
* Added new function for DEG volcano plot (plotDEGVolcano)
* Added new function for plotting pathway scores (plotPathway)
* Added Pathway Analysis section to the tutorials
* Added seed parameter to several functions and UI for reproducibility
* Updated R console and GUI tutorials to match each other
* Fixed console logging in the GUI

Changes in Version 2.4.0 (2021-10-27)
================================================================================
* Updated version to match Bioconductor

Changes in Version 2.3.2 (2021-10-24)
================================================================================

* Added summary table into the cellQC report
* Improved formatting in QC report
* Added functions getDEGTopTable() & plotBatchCorrCompare()
* Other refactors and bug fixes

Changes in Version 2.3.1 (2021-10-15)
================================================================================

* Several bug fixes

Changes in Version 2.2.2 (2021-10-10)
================================================================================

* Several enhancements, refactors, and bug fixes to the UI
* Refactor documentation and pkgdown site
* Added tutorials for R console analysis
* Updates to the UMAP generation in the SCTK-QC pipeline
* Addition of VAM to Pathway prediction tab
* Bug fix to the mitochondrial gene set functions

Changes in Version 2.1.3 (2021-05-14)
================================================================================

* Added diffAbundanceFET and plotClusterAbundance function
* Linked Shiny UI help buttons to new online help pages
* Several bug fixes

Changes in Version 2.0.2 (2021-05-08)
================================================================================

* Expanded convertSCEtoSeurat() function to copy additional data
* Updated and merged pkgdown docs
* Added HTML reports for Seurat curated workflow
* Refactor of Normalization UI
* Added generic wrapper function for dimensionality reduction
* Added tagging system for matrix type
* Several bug fixes
* Added missing documentation
* Added wrapper functions for normalization, dimensionality reduction and feature selection
* Added function seuratReport() to generate a seurat report from input SCE object

Changes in Version 2.0.1 (2021-01-07)
================================================================================

* Added cell type labeling functional, wrapping SingleR method
* Added cell type labeling UI under differential expression tab
* Added marker identification in Seurat workflow

Changes in Version 2.0.0 (2020-10-16)
================================================================================

* Added quality control (empty droplet detection, doublet detection, etc) functionality
* Ability to import data from varying preprocessing tools
* Ability to export SingleCellExperiment object as varying file types (flat file, Python anndata)
* Added functions for visualization of data
* New CellViewer functionality in UI
* Improvements to differential expression, now includes DESeq2, limma, ANOVA
* Incorporates Seurat workflow

Changes in Version 1.1.26 (2018-10-23)
================================================================================

* New UI design for the Differential Expression tab.
* New UI design for the Data Summary & Filtering tab.
* Support for additional assay modification including log transforming any assay and renaming assays.
* New function visPlot for creating scatterplots, boxplots, heatmaps, and barplots for custom gene sets.
* The Downsample tab now works on a generic counts matrix
* You can upload a SCtkExperiment object or a SingleCellExperiment object saved in an RDS file on the Upload tab.
* Differential Expression results can now be saved in the rowData of the object and loaded for later analysis.
* Improved ability to save a biomarker based on user options.
* The Differential Expression plot is not automatically created, for more user control with large datasets.

Changes in Version 1.1.3
================================================================================

* Improvements to plotting, change text size and hide labels in gsva plots.
* MAST violin and linear model plots are now more square when plotting less than 49 facets.
* Changed y axis label in plotBatchVariance to "Percent Explained Variation"

Changes in Version 1.1.2
================================================================================

* Ability to hide version number in the SCTK GUI.

Changes in Version 1.1.1
================================================================================

* Fixed a bug that would cause the diffex color bar to not display when special
characters were in the annotation.

Changes in Version 0.99.3
================================================================================

* Consistent use of camel case throughout package

Changes in Version 0.6.3
================================================================================

* Additional links to help documentation
* Example matrices on upload page.

Changes in Version 0.4.7
================================================================================

* Ability to download/reupload annotation data frame and convert annotations to
factors/numerics

Changes in Version 0.4.5
================================================================================

* Documentation updates to fix NOTES and pass BiocCheck
