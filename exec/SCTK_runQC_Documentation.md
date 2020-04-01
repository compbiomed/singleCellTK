# Generation of comprehensive quality control metrics with SCTK

This pipeline will import data from single-cell preprocessing algorithms (e.g. CellRanger), generate various quality control metrics (e.g. doublet scores), and output results in standard data containers (e.g. SingleCellExperiment).
Both the original droplet matrix and the filtered cell matrix will be processed.
This pipeline is focused on single cell data generated from microfluidic devices (e.g. 10X).

## Specifications

* The pipeline is currently written in the R language. Users will need to install R version 3.6.2 (or higher) in order to run all of the required packages. 
* For importing files from the HCA Optimus pipeline, the "scipy" module needs to be installed in the default version of Python on the system.

## Installation
The script will automatically try to install the "singleCellTK" package from Bioconductor if not available. However, currently this code is only located on a development branch which needs to be installed from GitHub:

```
library(devtools)
install_github("compbiomed/singleCellTK@importQC")
```

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

* [CellRanger (V2 or V3)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
* [STARsolo](https://github.com/alexdobin/STAR/releases)
* [Human Cell Atlas (HCA) Optimus pipeline](https://data.humancellatlas.org/pipelines/optimus-workflow)
* [BUStools](https://github.com/BUStools/bustools)
* [SEQC](https://github.com/ambrosejcarr/seqc)

For each tool, the pipeline expect the data in a certain format within a specific directory structure. For example, data generated with CellRanger V3 is expected to be under a sample folder (specified with the "-s" flag) and a base folder (specified with the "-b" flag).
Some tools prepend the sample name onto the output files instead of making a separate subdirectory (e.g. BUStools and SEQC). The combination of --base_path ("-b") and --sample ("-s") should specify the location of the data files:

![](/exec/SCTK_QC_Import.png)

## Arguments

The arguments are as follows:

-b, --base_path. Base path for the output from the preprocessing algorithm

-p, --preproc. Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"

-s, --sample. Name of the sample. This will be prepended to the cell barcodes.

-o, --directory. Output directory. A new subdirectory will be created with the name "sample". R, Python, and FlatFile directories will be created under the "sample" directory containing the data containers with QC metrics. Default ".".

-g, --gmt. GMT file containing gene sets for quality control. 

-t, --delim. Delimiter used in GMT file. Default "\t".

## Analyzing genes sets

Quantifying the level of gene sets can be useful quality control. For example, the percentage of counts from mitochondrial genes can be an indicator or cell stress or death. 

Users can pass a [GMT](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) file to the pipeline with one row for each gene set. The first column should be the name of the gene set (e.g. mito). 
The second column for each gene set in the GMT file (i.e. the description) should contain the location of where to look for the matching IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If a character string or an integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table.

## QC Outputs
This pipeline will run several QC algorithms. The QC metrics will be stored as `colData` within the outputed `singleCellExperiment` object in the R directory or in the file `colData.txt.gz` within the FlatFile directory. Here is a list of the currently available QC outputs:

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
| dropletUtils_barcodeRank_knee | Numeric, specifies total count at the knee point | dropletUtils |
| dropletUtils_barcodeRank_inflection | Numeric,  specifies total count at the inflection point | dropletUtils |

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
| decontX_contamination | Probability of contamination determined by decontX | celda |
| decontX_clusters | Clusters determined by Celda, a clustering algorithm that runs in the background of decontX | celda |


## Documentation of tools that are currently available within the pipeline:
#### Empty droplet detection:
* [emptyDrops](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html) from the package [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)

#### Doublet Detection
* [doubletCells](https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html) from the package [scran](http://bioconductor.org/packages/release/bioc/html/scran.html)
* [cxds](https://rdrr.io/bioc/scds/man/cxds.html), [bcds](https://rdrr.io/bioc/scds/man/bcds.html), and [cxds_bcds_hybrid](https://rdrr.io/bioc/scds/man/cxds_bcds_hybrid.html) from the package [scds](http://bioconductor.org/packages/release/bioc/html/scds.html)

#### Ambient RNA detection
* [decontX](https://rdrr.io/bioc/celda/man/decontX.html) from the package [celda](https://bioconductor.org/packages/release/bioc/html/celda.html)


