# Generation of comprehensive quality control metrics with SCTK

This pipeline will import data from single-cell preprocessing algorithms (e.g. CellRanger), generate various quality control metrics (e.g. doublet scores), and output results in standard data containers (e.g. SingleCellExperiment).
Both the original droplet matrix and the filtered cell matrix will be processed.
This pipeline is focused on single cell data generated from microfluidic devices (e.g. 10X).

## Specifications

* The pipeline is currently written in the R language. Users will need to install R version 3.6.2 (or higher) in order to run all of the required packages. The pipeline script "SCTK_runQC.R" can be download [here](https://github.com/rz2333/singleCellTK/blob/devel/exec/SCTK_runQC.R).
* For importing files from the HCA Optimus pipeline, the "scipy" module needs to be installed in the default version of Python on the system.
* The pipeline depends on some Python package for Python > 3.0.0. User will need to install Python 3.6.3 (or higher) in order to install the proper version of Python packages. 

## Installation
### Docker users
If you are using docker image to run the pipeline, please skip this section and refer sections **Running the pipeline** and **Docker and Singularity Images**. Docker image has all dependencies installed properly so that you can use it directly on any machine. 

### Install Python packages and dependencies
This pipeline depends on Python package described below. If these packages is not installed or not at the proper version, please install them by running the following code:

```
pip3 install --user h5py==2.9.0
pip3 install --user Scrublet
pip3 install --user virtualenv
pip3 install --user scanpy
```

### Install R packages and dependencies
The script will automatically try to install the "singleCellTK" package from Bioconductor if not available. However, currently this code is only located on a development branch which needs to be installed from GitHub:

```
library(devtools)
install_github("compbiomed/singleCellTK@devel")
```

## Data formats

This pipeline can currently import data from the following tools:

* [CellRanger (V2 or V3)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
* [STARsolo](https://github.com/alexdobin/STAR/releases)
* [Human Cell Atlas (HCA) Optimus pipeline](https://data.humancellatlas.org/pipelines/optimus-workflow)
* [BUStools](https://github.com/BUStools/bustools)
* [SEQC](https://github.com/ambrosejcarr/seqc)
* [SingleCellExperiment object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) saved in RDS file
* SingleCellExperiment object saved in [AnnData](https://github.com/theislab/anndata) hdf5 file
* Droplet / Cell matrix saved in txt/mtx files

For each tool, the pipeline expect the data in a certain format within a specific directory structure. For example, data generated with CellRanger V3 is expected to be under a sample folder (specified with the "-s" flag) and a base folder (specified with the "-b" flag).
Some tools prepend the sample name onto the output files instead of making a separate subdirectory (e.g. BUStools and SEQC). The combination of --base_path ("-b") and --sample ("-s") should specify the location of the data files:

![](/exec/SCTK_QC.png)

## Arguments

### Required arguments
The required arguments are as follows:

-b, --basePath (required). Base path for the output from the preprocessing algorithm

-P, --preproc (required). Algorithm used for preprocessing. One of One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus', 'DropEst', 'SceRDS', 'CountMatrix' and 'AnnData'. 

-s, --sample (required). Name of the sample. This will be prepended to the cell barcodes.

-o, --directory (required). Output directory. A new subdirectory will be created with the name "sample". R, Python, and FlatFile directories will be created under the "sample" directory containing the data containers with QC metrics. Default ".". More information about output directory structure is explained in "Understanding outputs" section below. 

-F, --outputFormat (required). The output format of this QC pipeline. Currently, it supports R (RDS), FlatFile, Python (AnnData) and HTAN (manifest files that meets HTAN requirement).

-S, --splitSample (required). Save a SingleCellExperiment object for each sample. Default is TRUE. If FALSE, the data of all samples will be combined into one SingleCellExperiment object and this object will be outputed.

### Optional arguments
The optional arguments are as follows. Their usage depend on type of data and user-defined behaviour. 

<details><summary>optional arguments</summary>
<p>
-g, --gmt (optional). GMT file containing gene sets for quality control. <br> <br>

-t, --delim (optional, required when -g is specified). Delimiter used in GMT file. Default "\t". <br> <br>

-G, --genome (optional). The name of genome reference. This is only required for CellRangerV2 data. <br> <br>

-y, --yaml (optional). YAML file used to specify parameters of QC functions called by singleCellTK QC pipeline. Please check "Specify parameters for QC algorithms" section for details. <br> <br>

-c, --cellData (optional). The full path of the RDS file or Matrix file of the cell matrix. This would be use only when --preproc is SceRDS or CountMatrix. <br> <br>

-r, --rawData (optional). The full path of the RDS file or Matrix file of the droplet matrix. This would be provided only when --preproc is SceRDS or CountMatrix. <br> <br>

-C, --cellPath (optional). The directory contains matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz files originally generated by 10x CellrangerV2 or CellrangerV3 (files in the filtered_feature_bc_matrix directory). This argument only works when --preproc is CellRangerV2 or CellRangerV3. Default is NULL. If 'base_path' is NULL, 'cellPath' or 'rawPath' should be specified. <br> <br>

-R, --rawPath (optional). The directory contains matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz files originally generated by 10x CellrangerV2 or CellrangerV3 (files in the raw_feature_bc_matrix directory). This argument only works when --preproc is CellRangerV2 or CellRangerV3. Default is NULL. If 'base_path' is NULL, 'cellPath' or 'rawPath' should be specified. <br> <br>

-d, --dataType. Type of data as input. Default is Both, which means taking both droplet and cell matrix as input. If set as 'Droplet', it will only processes droplet data. If set as 'Cell', it will only processes cell data. <br> <br>

-D, --detectCells. Detect cells from droplet matrix. Default is FALSE. This argument is only eavluated when -d is 'Droplet'. If set as TRUE, cells will be detected and cell matrixed will be subset from the droplet matrix. Also, quality control will be performed on the detected cell matrix. <br> <br>

-m, --cellDetectMethod. Methods to detect cells from droplet matrix. Default is 'EmptyDrops'. This argument is only eavluated when -D is 'TRUE'. Other options could be 'Knee' or 'Inflection'. More information is provided in the section below. <br> <br>

-n, --numCores. Number of cores used to run the pipeline. By default is 1. Parallel computing is enabled if -n is greater than 1. <br> <br>

-T, --parallelType. Type of parallel computing used for parallel computing. Parallel computing used in this pipeline depends on "BiocParallel" package. Default is 'MulticoreParam'. It can be 'MulticoreParam' or 'SnowParam'. This argument will be evaluated only when numCores > 1. <br> <br>
</p>
</details>

## Running the pipeline

### Methods to run pipeline on CellrangerV2 or CellrangerV3 data set
This pipeline enables different ways to import CellrangerV2/CellrangerV3 for flexibility. 

1. If the cellranger data set is saved in the default [cellranger output directory](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview), you can load the data by running following code: 

For CellRangerV3:
```
Rscript SCTK_runQC.R \
-b /basepath \
-P CellRangerV3 \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```

For CellRangerV2, the reference used by cellranger needs to be specified by -G/--genome:
```
Rscript SCTK_runQC.R \
-b /basepath \
-P CellRangerV2 \
-s SampleName \
-o Output_Directory \
-S TRUE \
-G hg19 \
-F R,Python,FlatFile,HTAN
```

As shown in the **Data formats**, -b specify the base path and usually it's the folder where you ran 10x cellranger count. -s specify the sample name, which has to be same as the name of the sample folder under the base folder. The folder layout would look like the following:

```
├── BasePath
└── SampleName 
    ├── outs
    |   ├── filtered_feature_bc_matrix
    |   |   ├── barcodes.tsv.gz
    |   |   ├── features.tsv.gz
    |   |   └── matrix.mtx.gz
    |   ├── raw_feature_bc_matrix
    |   |   ├── barcodes.tsv.gz
    |   |   ├── features.tsv.gz
    |   |   └── matrix.mtx.gz
        ...
```

2. If the cellranger count output have been moved out of the default [cellranger output directory](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview), you can specified input using -R and -C:
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

### Methods to run pipeline on data set stored in RDS file or matrix(txt or mtx)
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
-F R,Python,FlatFile,HTAN \
-r /path/to/matrix/file/droplet.txt \
-c path/to/matrix/file/cell.txt
```

### Methods to run pipeline on data set generated by other algorithms
If your data is preprocessed by other algorithms, you can load the data as shown in the figure above. Basically, the templated is shown below:

```
Rscript SCTK_runQC.R \
-b /base/path \
-P Preprocessing_Algorithm \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN
```

### Methods to run pipeline on only droplet count or cell count data
User can choose to run QC pipeline on only one of the droplet count or cell count matrix, instead of running on both. In this case, the pipeline will only take single input and perform quality control on it. One of the example is shown below:

```
Rscript SCTK_runQC.R \
-b /base/path \
-P Preprocessing_Algorithm \
-s SampleName \
-o Output_Directory \
-S TRUE \
-F R,Python,FlatFile,HTAN \
-d Droplet \
-D TRUE \
-m EmptyDrops
```

-d argument is used to specify which count matrix to used in the QC pipeline. Default is "Both", which means take both droplet and cell count data as input. 

If -d argument is set as "Droplet", the QC pipeline will only take droplet count matrix as input and perform quality control. You can choose whether to detect cells from the droplet matrix by setting -D as TRUE. If yes, cell count matrix will be detected and the pipeline will also perform quality control on this matrix and output the result. You could further define the method used to detect cells from droplet matrix by setting -m argument. -m could be one of "EmptyDrops", "Knee" or "Inflection". "EmptyDrops" will keep cells that pass the "runEmptyDrops" function test. "Knee" and "Inflection" will keep cells that pass the knee or inflection point returned from "runBarcodeRankDrops" function. 

If -d argument is set as "Cell", the QC pipeline will only take cell count matrix as input and perform quality control. A figure showing the analysis steps and outputs of different inputs is shown below:

![](/exec/Single_Input.png)

### Specify parameters for QC algorithms
User can specify parameters for QC algorithms in this pipeline with a yaml file (supplied with -y/--yamlFile argument). The current supported QC algorithms including doublet dection (bcds, cxds, cxds_bcds_hybrid, doubletFinder, doubletCells and scrublet), decontamination (decontX), emptyDrop detection (emptyDrops) and barcodeRankDrops (barcodeRanks). A summary of each function is shown below:
![](/exec/QC_yaml.png)

An example of QC parameters yaml file is shown below:
```
Params:   ### should not be omitted
  emptyDrops:           ### Name of the QC functions in the pipeline
    lower: 50           ### Parameter for this function
    niters: 5000        ### Parameter for this function
    testAmbient: True   ### Parameter for this function

  barcodeRanks: 
    lower: 50
```
The format of yaml file can be found [here](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html). The parameters should be consistent with the parameters of each QC function in singleCellTK package. Parameters that are not defined in this yaml file will use the default value.

### Analyzing genes sets

Quantifying the level of gene sets can be useful quality control. For example, the percentage of counts from mitochondrial genes can be an indicator or cell stress or death. 

Users can pass a [GMT](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) file to the pipeline with one row for each gene set. The first column should be the name of the gene set (e.g. mito). 
The second column for each gene set in the GMT file (i.e. the description) should contain the location of where to look for the matching IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If a character string or an integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table.

### Parallel computing

SCTK QC pipeline enables parallel computing to speed up the analysis. Parallel computing is enabled by setting -n/--numCores > 1. The -n/--numCores is used to set the number of cores used for the pipeline. 

The backend of parallel computing is supported by "BiocParallel" package. Therefore, users can select different types of parallel evaluation by setting -T/--parallelType argument. Default is "MulticoreParam". Currently, "MulticoreParam" and "SnowParam" is supported for -T argument. However, "MulticoreParam" is not supported by Windows system. Winows user can choose "SnowParam" as the backend of parallel computing. 


## Understanding outputs
The output directory is created under the path specified by -o/--directory argument. Each sample is stoed in the subdirectory (named by -s/--sample argument) within this output direcotry. Within each sample directory, the each output format will be separated into subdirectories. The output file hierarchy is shown below:
```
(root; output directory)
├── level3Meta.csv
├── level4Meta.csv
└── sample1 
    ├── R
    |   ├── sample1_Droplets.rds
    |   └── sample1_Cells.rds
    ├── Python
    |   ├── Droplets
    |   |   └── sample1.h5ad
    |   └── Cells
    |       └── sample1.h5ad
    ├── FlatFile
    |   ├── Droplets
    |   |   ├── assays
    |   |   |   └── sample1_counts.mtx.gz
    |   |   ├── metadata
    |   |   |   └── sample1_metadata.rds
    |   |   ├── sample1_colData.txt.gz
    |   |   └── sample1_rowData.txt.gz 
    |   └── Cells
    |       ├── assays
    |       |   └── sample1_counts.mtx.gz
    |       ├── metadata
    |       |   └── sample1_metadata.rds
    |       ├── reducedDims
    |       |   ├──sample1_decontX_UMAP.txt.gz
    |       |   ├──sample1_scrublet_TSNE.txt.gz
    |       |   └──sample1_scrublet_UMAP.txt.gz
    |       ├── sample1_colData.txt.gz
    |       └── sample1_rowData.txt.gz 
    └── sample1_QCparameters.yaml
```

## Docker and Singularity Images

singleCellTK is available to use with both Docker and Singularity. This container be used as a portal to the singleCellTK command line interface. The work was done in R. The singleCellTK docker image is available from [Docker Hub](https://hub.docker.com/r/campbio/sctk_qc).

### Docker

 If you have not used docker before, you can follow the instruction to install and set up docker in [Windows](https://docs.docker.com/docker-for-windows/), [Mac](https://docs.docker.com/docker-for-mac/) or [Linux](https://runnable.com/docker/install-docker-on-linux). 

The Docker image can be obtained by running: 
```
docker pull campbio/sctk_qc
```

Noted that the transcriptome data and GMT file needed to be accessible to the container via mounted volume. In the below example, mount volumn is enabled for accessing input and output directory using argument -v. To learn more about mounted volumes, please check out [this post](https://docs.docker.com/storage/volumes/). 

The usage of each argument is the same as running command line analysis. Here is an example code to perform quanlity control of CellRangerV3 data singleCellTK docker:

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
* [barcodeRanks](https://rdrr.io/github/MarioniLab/DropletUtils/man/barcodeRanks.html) from the package [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)

#### Doublet Detection
* [doubletCells](https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html) from the package [scran](http://bioconductor.org/packages/release/bioc/html/scran.html)
* [cxds](https://rdrr.io/bioc/scds/man/cxds.html), [bcds](https://rdrr.io/bioc/scds/man/bcds.html), and [cxds_bcds_hybrid](https://rdrr.io/bioc/scds/man/cxds_bcds_hybrid.html) from the package [scds](http://bioconductor.org/packages/release/bioc/html/scds.html)
* [doubletFinder](https://rdrr.io/github/chris-mcginnis-ucsf/DoubletFinder/man/doubletFinder.html) from the package [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
* [Scrublet](https://bioconda.github.io/recipes/scrublet/README.html) from the package [scrublet](https://github.com/allonkleinlab/scrublet)

#### Ambient RNA detection
* [decontX](https://rdrr.io/bioc/celda/man/decontX.html) from the package [celda](https://bioconductor.org/packages/release/bioc/html/celda.html)


