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
install_github("compbiomed/singleCellTK@devel")
```

## Running the pipeline

To run the pipeline script, users will need to download the 'SCTK_runQC.R' and run the following code:

```
Rscript SCTK_runQC.R \
-b /base/path \
-P Preprocessing_Algorithm \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```

or 

```
SCTK_runQC.R \
-b /base/path/ \
-P Preprocessing_Algorithm \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```

if the Rscript executable can be found in your environment with the command

```
/usr/bin/env Rscript --vanilla
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

### Required arguments
The required arguments are as follows:

-b, --base_path (required). Base path for the output from the preprocessing algorithm

-P, --preproc (required). Algorithm used for preprocessing. One of One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus', 'DropEst', 'SceRDS' and 'CountMatrix'. 

-s, --sample (required). Name of the sample. This will be prepended to the cell barcodes.

-o, --directory (required). Output directory. A new subdirectory will be created with the name "sample". R, Python, and FlatFile directories will be created under the "sample" directory containing the data containers with QC metrics. Default ".".

-F, --outputFormat (required). The output format of this QC pipeline. Currently, it supports R (RDS), FlatFile, Python (AnnData) and HTAN.

-S, --split_sample (required). Save SingleCellExperiment object for each sample. Default is TRUE. If FALSE, the data of all samples will be combined into one SingleCellExperiment object and this object will be outputed.

### Optional arguments
The optional arguments are as follows. Their usage depend on type of data and user-defined behaviour. 

-g, --gmt (optional). GMT file containing gene sets for quality control. 

-t, --delim (optional, required when -g is specified). Delimiter used in GMT file. Default "\t".

-r, --reference (optional). The name of genome reference. This is only required for CellRangerV2 data.

-G, --genome (optional). The name of genome reference. This is only required for CellRangerV2 data.

-y, --yaml (optional). YAML file containing parameters called by singleCellTK QC functions. Please check the following documentation for details. 

-c, --cell_data (optional). The full path of the RDS file or Matrix file of the cell matrix. This would be use only when --preproc is SceRDS or CountMatrix.

-r, --raw_data (optional). The full path of the RDS file or Matrix file of the droplet matrix. This would be provided only when --preproc is SceRDS or CountMatrix.

-C, --cell_data_path (optional). The directory contains cell matrix, which has gene symbol/ID as rownames and cell barcodes as colnames. Default is NULL. If 'base_path' is NULL, both 'cell_data_path' and 'raw_data_path' should be specified.

-R, --raw_data_path (optional). The directory contains droplet matrix, which has gene symbol/ID as rownames and droplet barcodes as colnames. Default is NULL. If 'base_path' is NULL, both 'cell_data_path' and 'raw_data_path' should also be specified.

## Methods to load CellrangerV2 or CellrangerV3 data set
This pipeline enables different ways to import CellrangerV2/CellrangerV3 for flexibility. 

1. If the cellranger data set is saved in the default [cellranger output directory](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview), you can load the data by running following code: 

For CellRangerV3:
```
Rscript SCTK_runQC.R \
-b /base/path \
-P CellRangerV3 \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```

For CellRangerV2, the reference used by cellranger needs to be specified by -G/--genome:
```
Rscript SCTK_runQC.R \
-b /base/path \
-P CellRangerV2 \
-s SampleName \
-o Output_Directory \
-S TRUE \
-G hg19 \
-F R,Python,FlatFile,HTAN
```

2. If the cell matrix and droplet matrix have been moved out of the cellranger output directory, you can specified directories using -R and -C:
```
Rscript SCTK_runQC.R \
-P CellRangerV2 \
-C /path/to/cell/matrix \
-R /path/to/droplet/matrix \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```
In this case, you must skip -b arguments and you can also skip -G argument for CellRangerV2 data.

## Methods to load data set stored in RDS file or txt file (matrix)
1. If you data in stored as a SingleCellExperiment in RDS file, singleCellTK also supports that type of input. To run quality control with RDS file as input, run the following code:
```
Rscript SCTK_runQC.R \
-P SceRDS \
-s Samplename \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN \
-r /path/to/rds/file/droplet.RDS \
-c /path/to/rds/file/cell.RDS
```

2. If your input is stored in txt file as a matrix, which has barcodes as colnames and genes as rownames, run the following code to start the quality control pipeline:
```
Rscript SCTK_runQC.R \
-P CountMatrix \
-s Samplename \
-o Output_Directory \
-S TRUE \
-r /path/to/matrix/file/droplet.txt \
-f path/to/matrix/file/cell.txt
```

## Analyzing genes sets

Quantifying the level of gene sets can be useful quality control. For example, the percentage of counts from mitochondrial genes can be an indicator or cell stress or death. 

Users can pass a [GMT](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) file to the pipeline with one row for each gene set. The first column should be the name of the gene set (e.g. mito). 
The second column for each gene set in the GMT file (i.e. the description) should contain the location of where to look for the matching IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If a character string or an integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table.

## Understanding outputs

## Docker and Singularity Images

singleCellTK is available to use with both Docker and Singularity. This container be used as a portal to the singleCellTK command line interface. The work was done in R. The singleCellTK docker image is available from [Docker Hub](https://hub.docker.com/r/campbio/sctk_qc).

### Docker

 If you have not used docker before, you can follow the instruction to install and set up docker in [Windows](https://docs.docker.com/docker-for-windows/), [Mac](https://docs.docker.com/docker-for-mac/) or [Linux](https://runnable.com/docker/install-docker-on-linux). 

The Docker image can be obtained by running: 
```
docker pull campbio/sctk_qc
```

Noted that the transcriptome data and GMT file needed to be accessible to the container via mounted volume. In the below example, mount is created for the input and output directory using argument -v. Here is an example code to perform quanlity control of CellRangerV3 data singleCellTK docker:

```
docker run --rm -v /path/to/data:/SCTK_docker \
-it campbio/sctk_qc:1.7.5 \
-b /SCTK_docker/cellranger \
-P CellRangerV3 \
-s pbmc_100x100 \
-o /SCTK_docker/result/tenx_v3_pbmc \
-g /SCTK_docker/mitochondrial_human_symbol.gmt \
-S TRUE \
-F R,Python,Flatfile,HTAN
```

The usage of each argument is the same as the original command line interface. To learn more about mounted volumes, please check out [this post](https://docs.docker.com/storage/volumes/).

### Singularity

The Singulatiry image can easily be built using Docker Hub as a source:

```
singularity pull docker://campbio/sctk_qc:1.7.5
```

The usage of singleCellTK Singularity image is very similar to that of Docker. In Singularity 3.0+, the mount volume is [automatically overlaid](https://singularity.lbl.gov/docs-mount). However, you can use argument --bind/-B to specify your own mount volume. The example is shown as below:

```
singularity run sctk_qc_1.7.5.sif \
-b ./cellranger \
-P CellRangerV3 \
-s pbmc_100x100 \
-o ./result/tenx_v3_pbmc \
-g ./mitochondrial_human_symbol.gmt
```

The code above assumed that the dataset is in your current directory, which is automatically mounted by Singularity. If you run Singularity image on BU SCC，it's recommended to re-set the home directory to mount. Otherwise, the container will load libraries in the SCC shared libraries, which might cause some conflicts. You can point to some "sanitized home" using argument [-H/--home](https://singularity.lbl.gov/faq#solution-1-specify-the-home-to-mount). Also, you might want to specify cpu architecture when run the docker on BU SCC using #$ -l cpu_arch=broadwell|haswell|skylake|cascadelake. Because the python packages are compiled by SIMD instructions that are only available on these two cpu architectures. 

## Documentation of tools that are currently available within the pipeline:
#### Empty droplet detection:
* [emptyDrops](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html) from the package [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)

#### Doublet Detection
* [doubletCells](https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html) from the package [scran](http://bioconductor.org/packages/release/bioc/html/scran.html)
* [cxds](https://rdrr.io/bioc/scds/man/cxds.html), [bcds](https://rdrr.io/bioc/scds/man/bcds.html), and [cxds_bcds_hybrid](https://rdrr.io/bioc/scds/man/cxds_bcds_hybrid.html) from the package [scds](http://bioconductor.org/packages/release/bioc/html/scds.html)
* [doubletFinder](https://rdrr.io/github/chris-mcginnis-ucsf/DoubletFinder/man/doubletFinder.html) from the package [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
* [Scrublet](https://bioconda.github.io/recipes/scrublet/README.html) from the package [scrublet](https://github.com/allonkleinlab/scrublet)

#### Ambient RNA detection
* [decontX](https://rdrr.io/bioc/celda/man/decontX.html) from the package [celda](https://bioconductor.org/packages/release/bioc/html/celda.html)


