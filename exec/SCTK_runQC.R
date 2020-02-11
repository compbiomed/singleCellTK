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
option_list <- list(optparse::make_option(c("-b", "--base_path"),
        type="character",
        default=NULL,
        help="Base path for the output from the preprocessing algorithm"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRangerV3",
        help="Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"),
    optparse::make_option(c("-s","--sample"),
        type="character",
        help="Name of the sample. This will be prepended to the cell barcodes."),
    optparse::make_option(c("-o","--directory"),
        type="character",
        default=".",
        help="Output directory"),
    optparse::make_option(c("-g","--gmt"),
        type="character",
        default=NULL,
        help="GMT file containing gene sets for quality control. The second column in the GMT file (i.e. the description) should contain the location to look for the IDs in the data. If set to 'rownames', then the gene set IDs will be matched with the row IDs of the data matrix. If another character or integer index is supplied, then gene set IDs will be matched to IDs the that column of feature table."),
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

## Use appropriate import function for preprocessing tool
dropletSCE <- NULL
filteredSCE <- NULL
if (preproc == "BUStools") {
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
}

## Read in gene sets for QC
geneSetCollection <- NULL
if(!is.null(gmt)) {
  geneSetCollection <- GSEABase::getGmt(gmt, sep=sep)
}

## Run QC functions
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
  filteredSCE <- runCellQC(sce = filteredSCE, geneSetCollection = geneSetCollection)
}  


## Merge singleCellExperiment objects
mergedDropletSCE <- NULL
mergedFilteredSCE <- NULL
if(!is.null(filteredSCE) & !is.null(dropletSCE)) {
  mergedDropletSCE <- mergeSCEColData(dropletSCE, filteredSCE)
  mergedFilteredSCE <- mergeSCEColData(filteredSCE, dropletSCE)
} else {
  mergedFilteredSCE <- filteredSCE
}

## Create directories and save objects
dir.create(file.path(directory, samplename), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "R"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "Python"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "FlatFile"), showWarnings = TRUE, recursive = TRUE)

if(!is.null(mergedDropletSCE)){
  ## Export to R 
  fn <- file.path(directory, samplename, "R", paste0(samplename , "_Droplets.rds"))
  saveRDS(object = mergedDropletSCE, file = fn)
  
  ## Export to flatfile
  fn <- file.path(directory, samplename, "FlatFile", "Droplets")
  writeSCE(mergedDropletSCE, outputDir = fn)
}
if(!is.null(mergedFilteredSCE)) {
  ## Export to R    
  fn <- file.path(directory, samplename, "R", paste0(samplename , "_FilteredCells.rds"))
  saveRDS(object = mergedFilteredSCE, file = fn)

  ## Export to flatfile  
  fn <- file.path(directory, samplename, "FlatFile", "FilteredCells")
  writeSCE(mergedFilteredSCE, outputDir = fn)
}  

## ToDo ##
## Export to Python



sessionInfo()
