# Generation of comprehensive quality control metrics with SCTK

This pipeline will import data from single-cell preprocessing algorithms (e.g. CellRanger), generate various quality control metrics (e.g. doublet scores), and output results in standard data containers (e.g. SingleCellExperiment).
Both the original droplet matrix and the filtered cell matrix will be processed.
This pipeline is focused on single cell data generated from microfluidic devices (e.g. 10X).

## Specifications

* The pipeline is currently written in the R language. Users will need to install R version 3.6.2 (or higher) in order to run all of the required packages. 
* For importing files from the HCA Optimus pipeline, the "scipy" module needs to be installed in the default version of Python on the system.


## Running the pipeline

To run the pipeline script, users will need to download the 'SCTK_runQC.R' and run the following code:

```
Rscript SCTK_runQC.R -p /base/path -p Preprocessing_Algorithm -s SampleName -o Output_Directory
```

or 

```
SCTK_runQC.R -p /base/path/ -p Preprocessing_Algorithm -s SampleName -o Output_Directory
```

if the Rscript executable can be found in your environment with the command

```
/usr/bin/env Rscript
```


## Data formats

This pipeline can currently import data from the following tools:

* CellRanger (V2 or V3)
* STARsolo 
* Human Cell Atlas (HCA) Optimus pipeline
* BUStools
* SEQC

For each tool, the pipeline expect the data in a certain format within a specific directory structure. For example, data generated with CellRanger V3 is expected to be under a sample folder (specified with the "-s" flag) and a base folder (specified with the "-b" flag).
Some tools prepend the sample name onto the output files instead of making a separate subdirectory (e.g. BUStools and SEQC). The combination of --base_path ("-b") and --sample ("-s") should specify the location of the data files:

![](/exec/SCTK_QC_Import.png)

## Arguments

The arguments are as follows:

-b, --base_path. Base path for the output from the preprocessing algorithm

-p, --preproc. Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"

-s, --sample. Name of the sample. This will be prepended to the cell barcodes.

-o, --directory. Output directory. A new subdirectory will be created with the name "sample". R, Python, and FlatFile directories will be created under the "sample" directory containing the data containers with QC metrics. Default ".".

-g, --gmt. GMT file containing gene sets for quality control. The second column in the GMT file (i.e. the description) should contain the location to look for the IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If another character or integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table.

-t, --delim. Delimiter used in GMT file. Default "\t".

## Including genes sets for analysis

## Documentation of tools that are currently available within the pipeline:
#### Empty droplet detection:
- emptyDrops(dropletUtils): https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html

#### Doublet Detection:
- doubletCells(scran): https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html


