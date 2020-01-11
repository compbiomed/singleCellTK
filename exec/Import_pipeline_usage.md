# importScript documentation

This pipeline will import outputs from single-cell preprocessing steps (i.e. CellRanger, STARSolo, BUStools), run various quality control metrics (ex. doublet detection), and output SingleCellExperiment object(s) containing QC metric information.

## Specifications

The pipeline is currently written in the R language. Users will need to install R version 3.6.0 (or higher) in order to run all of the required software. 

## Running the import pipeline

To run the pipeline script, users will need to upload the importScript.R Rscript to the desired folder and run the following code:

```
Rscript singleCellTK_SampleQC.R -u /path/to/raw_feature_bc_matrix -f /path/to/filtered_feature_bc_matrix -p Preprocessing -g TRUE -s SampleName -d Directory
```

The arguments are as follows:

-d The path to the unfiltered/raw output from preprocessing steps. A "matrix.mtx" file containing the counts data, "features.tsv"containing the features data, and a "barcodes.tsv" containing the barcodes for all of the samples is required.

-c The path to the filtered output from preprocessing steps. A "matrix.mtx" file containing the counts data, "features.tsv"containing the features data, and a "barcodes.tsv" containing the barcodes for all of the samples is required.

-p Preprocessing step used (CellRanger, etc.)

-g Whether the outputs from -u/-f are gzipped or not

-s The desired sample name

-d The desired output directory name

If users wish to run this script as a batch job on a computing cluster, we have provided a .sh script that may be used to run the pipeline Rscript:

```
qsub singleCellTK_SampleQC.sh
```

## Documentation of tools that are currently available within the pipeline:
#### Empty droplet detection:
- emptyDrops(dropletUtils): https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html

#### Doublet Detection:
- doubletCells(scran): https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html


