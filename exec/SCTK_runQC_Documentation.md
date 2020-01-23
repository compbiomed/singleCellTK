# SCTK QC pipeline

This pipeline will import outputs from single-cell preprocessing algorithms (e.g. CellRanger), generate various quality control metrics (e.g. doublet scores), and output SingleCellExperiment objects containing QC metric information.

## Specifications

The pipeline is currently written in the R language. Users will need to install R version 3.6.2 (or higher) in order to run all of the required software. 
For importing files from the HCA Optimus pipeline, the "scipy" module needs to be installed in the default version of Python.


## Running the pipeline

To run the pipeline script, users will need to upload the importScript.R Rscript to the desired folder and run the following code:

```
Rscript SCTK_runQC.R -p /path/to/data -p Preprocessing_Algorithm -s SampleName -o Output_Directory
```

or 

```
SCTK_runQC.R -p /path/to/data -p Preprocessing_Algorithm -s SampleName -o Output_Directory
```

if Rscript can be found in your environment with the command

```
/usr/bin/env Rscript
```


## Data formats

![](exec/SCTK_QC_Import.png)

## Arguments

The arguments are as follows:

-d The path to the unfiltered/raw output from preprocessing steps. A "matrix.mtx" file containing the counts data, "features.tsv"containing the features data, and a "barcodes.tsv" containing the barcodes for all of the samples is required.

-c The path to the filtered output from preprocessing steps. A "matrix.mtx" file containing the counts data, "features.tsv"containing the features data, and a "barcodes.tsv" containing the barcodes for all of the samples is required.

-p Preprocessing step used (CellRanger, etc.)

-g Whether the outputs from -u/-f are gzipped or not

-s The desired sample name

-o The desired output directory name

## Documentation of tools that are currently available within the pipeline:
#### Empty droplet detection:
- emptyDrops(dropletUtils): https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html

#### Doublet Detection:
- doubletCells(scran): https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html


