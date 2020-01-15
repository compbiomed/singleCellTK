#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("singleCellTK"))

##Check to see if necessary packages are installed

#Bioconductor packages
bioc.packages <- c("singleCellTK")

bioc.package.check <- lapply(bioc.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        BiocManager::install("singleCellTK")
    }
})

#CRAN packages
cran.packages <- c()

cran.package.check <- lapply(cran.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
    }
})


##Read in flags from command line using optparse

option_list <- list(optparse::make_option(c("-d", "--droplet"),
        type="character",
        default=NA,
        help="path to the unfiltered output from preprocessing steps [CellRanger, etc.]"),
    optparse::make_option(c("-c", "--cell"),
        type="character",
        default=NA,
        help="path to the filtered output from preprocessing steps [BUStools, etc.]"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRanger",
        help="One of 'CellRanger', 'BUStools', 'STARSolo'"),
    optparse::make_option(c("-g","--gzip"),
        type="logical",
        default=TRUE,
        help="Are your matrix, barcode, and features files gzipped?"),
    optparse::make_option(c("-s","--samplename"),
        type="character",
        help="Sample name"),
    optparse::make_option(c("-o","--directory"),
        type="character",
        default="R",
        help="Directory for output SingleCellExperiment objects"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
droplet.path <- opt$droplet
filtered.path <- opt$cell
preproc <- opt$preproc
gzip <- opt$gzip
samplename <- opt$samplename
directory <- opt$directory

if (is.na(samplename)){
  stop("A sample name is required. Please specify using the -s flag.")
}

## Use appropriate import function for preprocessing tool

if (preproc == "BUStools") {
    dropletSCE <- importBUStools(BUStoolsDir = droplet.path, sample = "", gzipped = gzip)
    if(!is.na(filtered.path)){
      filteredSCE <- importBUStools(BUStoolsDir = filtered.path, sample = "", gzipped = gzip)
    }
} else if(preproc == "STARSolo"){
    dropletSCE <- importSTARsolo(STARsoloDir = droplet.path, sample = "", STARsoloOuts = "", gzipped = gzip, class = "Matrix")
    if(!is.na(filtered.path)){
      filteredSCE <- importSTARsolo(STARsoloDir = filtered.path, sample = "", STARsoloOuts = "", gzipped = gzip)
    }
} else if(preproc == "CellRanger"){
    dropletSCE <- importCellRanger(cellRangerDirs = droplet.path, samples = "", cellRangerOuts = "", gzipped = gzip, class = "Matrix")
    if(!is.na(filtered.path)){
      filteredSCE <- importCellRanger(cellRangerDirs = filtered.path, samples = "", cellRangerOuts = "", gzipped = gzip)
    }
} else {
  stop(paste0("'", preproc, "' not supported."))
}



## Run QC functions
dropletSCE <- runDropletQC(sce = dropletSCE)

if(!is.na(filtered.path)){
  filteredSCE <- runCellQC(sce = filteredSCE)
}

## Merge singleCellExperiment objects
if(!is.na(filtered.path)){
  mergedDropletSCE <- mergeSCEColData(dropletSCE, filteredSCE)
  mergedFilteredSCE <- mergeSCEColData(filteredSCE, dropletSCE)
}else{
  mergedDropletSCE <- dropletSCE[,!is.na(dropletSCE$dropletUtils_emptyDrops_fdr)]
  mergedDropletSCE <- mergedDropletSCE[,mergedDropletSCE$dropletUtils_emptyDrops_fdr < 0.25]
}

## Create directory
if(is.null(directory)){
  directory <- samplename
}
dir.create(file.path(directory), showWarnings = TRUE)

setwd(file.path(directory))

## Save singleCellExperiment object
dir.create(file.path("R"), showWarnings = TRUE)
setwd("R")
saveRDS(object = mergedDropletSCE, file = paste0(samplename , "_Droplets.rds"))

if(!is.na(filtered.path)){
  saveRDS(object = mergedFilteredSCE, file = paste0(samplename , "_FilteredCells.rds"))
}

## ToDo ##
## Export to Python
## Export to flatfile


sessionInfo()
