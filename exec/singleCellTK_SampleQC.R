#!/usr/bin/env Rscript

##########


##########


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
        help="Path to the unfiltered output from preprocessing steps, including all droplets. [CellRanger, etc.]"),
    optparse::make_option(c("-c", "--cell"),
        type="character",
        default=NA,
        help="Path to the filtered output from preprocessing steps. [CellRanger, etc.] Optional."),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRanger",
        help="One of 'CellRanger', 'BUStools', or 'STARSolo"),
    optparse::make_option(c("-g","--gzip"),
        type="logical",
        default=TRUE,
        help="Are your matrix, barcode, and features files gzipped?"),
    optparse::make_option(c("-s","--samplename"),
        type="character",
        help="Sample name"),
    optparse::make_option(c("-o","--directory"),
        type="character",
        default=NULL,
        help="Output directory for SingleCellExperiment objects created via the pipeline."))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

droplet.path <- opt$droplet

cell.path <- opt$cell

preproc <- opt$preproc

gzip <- opt$gzip

samplename <- opt$samplename

directory <- opt$directory

##Use appropriate import function for preprocessing tool
if(preproc == "BUStools") {
    dropletSCE <- importBUStools(BUStoolsDir=droplet.path,sample="",gzipped=gzip)
    cellSCE <- importBUStools(BUStoolsDir=cell.path,sample="",gzipped=gzip)
}else if(preproc == "STARSolo"){
    dropletSCE <- importSTARsolo(STARsoloDir=droplet.path,sample="",STARsoloOuts="",gzipped=gzip)
    cellSCE <- importSTARsolo(STARsoloDir=cell.path,sample="",STARsoloOuts="",gzipped=gzip)
}else if(preproc == "CellRanger"){
    dropletSCE <- importCellRanger(cellRangerDirs=droplet.path,samples="", cellRangerOuts="", gzipped=gzip)
    cellSCE <- importCellRanger(cellRangerDirs=cell.path,samples="", cellRangerOuts="", gzipped=gzip, class =)}

##Run Appropriate QC functions
dropletSCE = runQC(sce = dropletSCE, algorithms = "emptyDrops")
cellSCE = runQC(sce = cellSCE, algorithms = "doubletCells")

#Merge singleCellExperiment objects
mergedDropletSCE <- mergeSCEColData(dropletSCE, cellSCE)
mergedCellSCE <- mergeSCEColData(cellSCE, dropletSCE)

#Create directory
if(is.null(directory)){
    directory <- samplename
}

dir.create(file.path(directory), showWarnings = TRUE)

setwd(file.path(directory))

#Save singleCellExperiment object
dir.create(file.path("R"), showWarnings = TRUE)
setwd("R")
saveRDS(object = mergedDropletSCE, file = paste0(samplename , "_Droplets.rds"))
saveRDS(object = mergedCellSCE, file = paste0(samplename , "_FilteredCells.rds"))

sessionInfo()
