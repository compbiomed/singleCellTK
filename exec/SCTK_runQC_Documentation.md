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

-P, --preproc. Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"

-s, --sample. Name of the sample. This will be prepended to the cell barcodes.

-o, --directory. Output directory. A new subdirectory will be created with the name "sample". R, Python, and FlatFile directories will be created under the "sample" directory containing the data containers with QC metrics. Default ".".

-g, --gmt. GMT file containing gene sets for quality control. 

-t, --delim. Delimiter used in GMT file. Default "\t".

-r, --reference. The name of genome reference. This is only required for CellRangerV2 data.

## Analyzing genes sets

Quantifying the level of gene sets can be useful quality control. For example, the percentage of counts from mitochondrial genes can be an indicator or cell stress or death. 

Users can pass a [GMT](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) file to the pipeline with one row for each gene set. The first column should be the name of the gene set (e.g. mito). 
The second column for each gene set in the GMT file (i.e. the description) should contain the location of where to look for the matching IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If a character string or an integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table.

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
singularity run sctk_qc_0.1.5.sif \
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


