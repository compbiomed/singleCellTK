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

option_list <- list(optparse::make_option(c("-u", "--unfiltered"),
        type="character",
        default=NA,
        help="path to the unfiltered output from preprocessing steps [CellRanger, etc.]"),
    optparse::make_option(c("-f", "--filtered"),
        type="character",
        default=NA,
        help="path to the filtered output from preprocessing steps [BUStools, etc.]"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRanger",
        help="e.g. CellRanger, BUStools, STARSolo"),
    optparse::make_option(c("-g","--gzip"),
        type="logical",
        default=TRUE,
        help="Are your matrix, barcode, and features files gzipped?"),
    optparse::make_option(c("-s","--samplename"),
        type="character",
        default=NA,
        help="Sample name"),
    optparse::make_option(c("-d","--directory"),
        type="character",
        default="R",
        help="Directory for output SingleCellExperiment objects"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

unfiltered.path <- opt$unfiltered

filtered.path <- opt$filtered

preproc <- opt$preproc

gzip <- opt$gzip

samplename <- opt$samplename

directory <- opt$directory

##Use appropriate import function for preprocessing tool
if(preproc == "BUStools") {
    unfilteredSCE <- importBUStools(BUStoolsDir=unfiltered.path,sample="",gzipped=gzip, class="Matrix")
    filteredSCE <- importBUStools(BUStoolsDir=filtered.path,sample="",gzipped=gzip, class="Matrix")
}else if(preproc == "STARSolo"){
    unfilteredSCE <- importSTARsolo(STARsoloDir=unfiltered.path,sample="",STARsoloOuts="",gzipped=gzip, class="Matrix")
    filteredSCE <- importSTARsolo(STARsoloDir=filtered.path,sample="",STARsoloOuts="",gzipped=gzip, class="Matrix")
}else if(preproc == "CellRanger"){
    unfilteredSCE <- importCellRanger(cellRangerDirs=unfiltered.path,samples="", cellRangerOuts="", gzipped=gzip, class="Matrix")
    filteredSCE <- importCellRanger(cellRangerDirs=filtered.path,samples="", cellRangerOuts="", gzipped=gzip, class ="Matrix")}

##Run Appropriate QC functions
unfilteredSCE = runQC(sce = unfilteredSCE, algorithms = "emptyDrops")
filteredSCE = runQC(sce = filteredSCE, algorithms = "doubletCells")

#Merge singleCellExperiment objects
mergedUnfilteredSCE <- mergeSCEColData(unfilteredSCE, filteredSCE)
mergedFilteredSCE <- mergeSCEColData(filteredSCE, unfilteredSCE)

#Create directory
dir.create(file.path(directory), showWarnings = TRUE)
setwd(file.path(directory))

#Save singleCellExperiment object
saveRDS(object = mergedUnfilteredSCE, file = paste0(samplename , "_Droplets.rds"))
saveRDS(object = mergedFilteredSCE, file = paste0(samplename , "_FilteredCells.rds"))

sessionInfo()
