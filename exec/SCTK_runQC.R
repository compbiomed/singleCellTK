#!/usr/bin/env Rscript

##Check to see if necessary packages are installed
#CRAN packages
cran.packages <- c("optparse")

cran.package.check <- lapply(cran.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
    }
})

#Bioconductor packages
bioc.packages <- c("singleCellTK")

bioc.package.check <- lapply(bioc.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        BiocManager::install("singleCellTK")
    }
})


##Read in flags from command line using optparse
<<<<<<< HEAD
require("optparse")
require("singleCellTK")
option_list <- list(optparse::make_option(c("-d", "--droplet"),
        type="character",
        default=NA,
        help="path to the unfiltered droplet counts matrix"),
    optparse::make_option(c("-c", "--cell"),
        type="character",
        default=NA,
        help="path to the filtered cells counts matrix"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRanger",
        help="One of 'CellRanger', 'BUStools', 'STARSolo', 'SEQC', or 'Optimus'"),
    optparse::make_option(c("-g","--gzip"),
        type="logical",
        default=TRUE,
        help="Are your matrix, barcode, and features files gzipped?"),
    optparse::make_option(c("-s","--samplename"),
=======
option_list <- list(optparse::make_option(c("-b", "--base_path"),
        type="character",
        default=NULL,
        help="Base path for the output from the preprocessing algorithm"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRangerV3",
        help="Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"),
    optparse::make_option(c("-s","--sample"),
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
        type="character",
        help="Name of the sample. This will be prepended to the cell barcodes."),
    optparse::make_option(c("-o","--directory"),
        type="character",
        default=".",
        help="Output directory"),
    optparse::make_option(c("-g","--gmt"),
        type="character",
        default=NULL,
<<<<<<< HEAD
        help="GMT file containing gene sets for quality control"),
=======
        help="GMT file containing gene sets for quality control. The second column in the GMT file (i.e. the description) should contain the location to look for the IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If another character or integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table."),
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
    optparse::make_option(c("-t","--delim"),
        type="character",
        default="\t",
        help="Delimiter used in GMT file"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
path <- opt$base_path
preproc <- opt$preproc
samplename <- opt$sample
directory <- opt$directory
gmt <- opt$gmt
sep <- opt$delim

<<<<<<< HEAD
if (is.na(samplename)){
  stop("A sample name is required. Please specify using the -s flag.")
}

if(is.na(droplet.path) && is.na(filtered.path)){
  stop("Either the droplet counts or the filtered counts file path need to be specified.")
}

=======
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
## Use appropriate import function for preprocessing tool
dropletSCE <- NULL
filteredSCE <- NULL
if (preproc == "BUStools") {
<<<<<<< HEAD
    if(!is.na(droplet.path)){
        dropletSCE <- importBUStools(BUStoolsDir = droplet.path, sample = "", gzipped = gzip)
    }
    if(!is.na(filtered.path)){
        filteredSCE <- importBUStools(BUStoolsDir = filtered.path, sample = "", gzipped = gzip)
    }
} else if(preproc == "STARSolo"){
    if(!is.na(droplet.path)){
        dropletSCE <- importSTARsolo(STARsoloDir = droplet.path, sample = "", STARsoloOuts = "", gzipped = gzip, class = "Matrix")
    }
    if(!is.na(filtered.path)){
        filteredSCE <- importSTARsolo(STARsoloDir = filtered.path, sample = "", STARsoloOuts = "", gzipped = gzip)
    }
} else if(preproc == "CellRanger"){
    if(!is.na(droplet.path)){
        dropletSCE <- importCellRanger(cellRangerDirs = droplet.path, samples = "", cellRangerOuts = "", gzipped = gzip, class = "Matrix")
    }
    if(!is.na(filtered.path)){
        filteredSCE <- importCellRanger(cellRangerDirs = filtered.path, samples = "", cellRangerOuts = "", gzipped = gzip)
    }
} else if(preproc == "SEQC"){
    if(!is.na(droplet.path)){
        dropletSCE <- importSEQC(seqcDirs = droplet.path, samples = samplename, prefix = samplename, gzipped = gzip, class = "Matrix")
    }
    if(!is.na(filtered.path)){
        filteredSCE <- importSEQC(seqcDirs = filtered.path, samples = samplename, prefix = samplename, gzipped = gzip)
    }
} else if(preproc == "Optimus"){
    if(!is.na(droplet.path)){
        dropletSCE <- importOptimus(OptimusDirs = droplet.path, samples = samplename)
    }
    if(!is.na(filtered.path)){
        filteredSCE <- importOptimus(OptimusDirs = droplet.path, samples = samplename)
    }
} else{
    stop(paste0("'", preproc, "' not supported."))
=======
  dropletSCE <- importBUStools(BUStoolsDir = path, sample = samplename, class = "Matrix")
} else if(preproc == "STARSolo"){
  dropletSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/raw", class = "Matrix")
  filteredSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix")
} else if(preproc == "CellRangerV3"){
  dropletSCE <- importCellRanger(cellRangerDirs = path, samples = samplename, gzipped = TRUE, cellRangerOuts = "outs/raw_feature_bc_matrix", class = "Matrix")
  filteredSCE <- importCellRanger(cellRangerDirs = path, samples = samplename, gzipped = TRUE, cellRangerOuts = "outs/filtered_feature_bc_matrix", class = "Matrix")
} else if(preproc == "CellRangerV2"){
  dropletSCE <- importCellRanger(cellRangerDirs = path, samples = samplename, gzipped = FALSE, cellRangerOuts = "outs/raw_feature_bc_matrix", class = "Matrix")
  filteredSCE <- importCellRanger(cellRangerDirs = path, samples = samplename, gzipped = FALSE, cellRangerOuts = "outs/filtered_gene_bc_matrix", class = "Matrix")
} else if(preproc == "SEQC"){
  dropletSCE <- importSEQC(seqcDirs = path, samples = samplename, prefix = samplename, class = "Matrix")
} else if(preproc == "Optimus"){
  dropletSCE <- importOptimus(OptimusDirs = path, samples = samplename)
  filteredSCE <- dropletSCE[,which(dropletSCE$dropletUtils_emptyDrops_IsCell)]
} else {
  stop(paste0("'", preproc, "' not supported."))
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
}

## Read in gene sets for QC
geneSetCollection <- NULL
if(!is.null(gmt)) {
  geneSetCollection <- GSEABase::getGmt(gmt, sep=sep)
}

## Run QC functions
<<<<<<< HEAD
if(!is.na(droplet.path)){
  dropletSCE <- runDropletQC(sce = dropletSCE)
}
if(!is.na(filtered.path)){
=======
if(!is.null(dropletSCE)) {
  message(paste0(date(), " .. Running droplet QC"))    
  dropletSCE <- runDropletQC(sce = dropletSCE)
  
  if(is.null(filteredSCE)) {
    ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
    filteredSCE <- dropletSCE[,ix]
  }  
}

if(!is.null(filteredSCE)) {
  message(paste0(date(), " .. Running cell QC"))    
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
  filteredSCE <- runCellQC(sce = filteredSCE, geneSetCollection = geneSetCollection)
}  


## Merge singleCellExperiment objects
<<<<<<< HEAD
if(!is.na(filtered.path) && !is.na(droplet.path)){
  mergedDropletSCE <- mergeSCEColData(dropletSCE, filteredSCE)
  mergedFilteredSCE <- mergeSCEColData(filteredSCE, dropletSCE)
}else if(is.na(filtered.path)){
  mergedDropletSCE <- dropletSCE[,!is.na(dropletSCE$dropletUtils_emptyDrops_fdr)]
  mergedDropletSCE <- mergedDropletSCE[,mergedDropletSCE$dropletUtils_emptyDrops_fdr < 0.05]
}else if(is.na(droplet.path)){
  mergedFilteredSCE <- filteredSCE
}

## Create directory
if(is.null(directory)){
  directory <- samplename
=======
mergedDropletSCE <- NULL
mergedFilteredSCE <- NULL
if(!is.null(filteredSCE) & !is.null(dropletSCE)) {
  mergedDropletSCE <- mergeSCEColData(dropletSCE, filteredSCE)
  mergedFilteredSCE <- mergeSCEColData(filteredSCE, dropletSCE)
} else {
  mergedFilteredSCE <- filteredSCE
>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b
}

## Create directories and save objects
dir.create(file.path(directory, samplename), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "R"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "Python"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "FlatFile"), showWarnings = TRUE, recursive = TRUE)

<<<<<<< HEAD
## Save singleCellExperiment object
dir.create(file.path("R"), showWarnings = TRUE)
setwd("R")

if(!is.na(filtered.path)){
  saveRDS(object = mergedFilteredSCE, file = paste0(samplename , "_FilteredCells.rds"))
}
if(!is.na(droplet.path)){
  saveRDS(object = mergedDropletSCE, file = paste0(samplename , "_Droplets.rds"))
}
=======
if(!is.null(mergedDropletSCE)){
  fn <- file.path(directory, samplename, "R", paste0(samplename , "_Droplets.rds"))
  saveRDS(object = mergedDropletSCE, file = fn)
}
if(!is.null(mergedFilteredSCE)) {
  fn <- file.path(directory, samplename, "R", paste0(samplename , "_FilteredCells.rds"))
  saveRDS(object = mergedFilteredSCE, file = fn)
}  

>>>>>>> 58756e23cbb48bc90cc1ee98a17815a03230bc2b

## ToDo ##
## Export to Python
## Export to flatfile


sessionInfo()
