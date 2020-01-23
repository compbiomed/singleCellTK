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
option_list <- list(optparse::make_option(c("-d", "--droplet"),
        type="character",
        default=NULL,
        help="Base path for the unfiltered droplet counts matrix"),
    optparse::make_option(c("-c", "--cell"),
        type="character",
        default=NULL,        
        help="Base path for the filtered cells counts matrix"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRanger",
        help="Algorithm used for preprocessing. One of 'CellRanger', 'BUStools', 'STARSolo', 'SEQC', 'Optimus'"),
    optparse::make_option(c("-s","--sample"),
        type="character",
        help="Sample name"),
    optparse::make_option(c("-o","--directory"),
        type="character",
        default=".",
        help="Output directory"),
    optparse::make_option(c("-g","--gmt"),
        type="character",
        default=NULL,
        help="GMT file containing gene sets for quality control"),
    optparse::make_option(c("-t","--delim"),
        type="character",
        default="\t",
        help="Delimiter used in GMT file"))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
droplet.path <- opt$droplet
filtered.path <- opt$cell
preproc <- opt$preproc
samplename <- opt$sample
directory <- opt$directory
gmt <- opt$gmt
sep <- opt$delim

## Check parameters
if (is.null(samplename)){
  stop("A sample name is required. Please specify using the -s flag.")
}

if(is.null(droplet.path) && is.null(filtered.path)){
  stop("Either the droplet counts or the filtered counts file path need to be specified.")
}

## Use appropriate import function for preprocessing tool
if (preproc == "BUStools") {
    if(!is.null(droplet.path)){
      dropletSCE <- importBUStools(BUStoolsDir = droplet.path, sample = samplename, class = "Matrix")
    }  
    if(!is.null(filtered.path)){
      filteredSCE <- importBUStools(BUStoolsDir = filtered.path, sample = samplename, class = "Matrix")
    }
} else if(preproc == "STARSolo"){
    if(!is.null(droplet.path)){
      dropletSCE <- importSTARsolo(STARsoloDir = droplet.path, sample = samplename, STARsoloOuts = "/outs/raw_feature_bc_matrix", class = "Matrix")
      dropletSCE$sample <- samplename
    }  
    if(!is.null(filtered.path)){
      filteredSCE <- importSTARsolo(STARsoloDir = filtered.path, sample = samplename, STARsoloOuts = "/outs/filtered_feature_bc_matrix", class = "Matrix")
    }
} else if(preproc == "CellRanger"){
    if(!is.null(droplet.path)){
      dropletSCE <- importCellRanger(cellRangerDirs = droplet.path, samples = samplename, cellRangerOuts = "/outs/raw_feature_bc_matrix", class = "Matrix")
    }  
    if(!is.null(filtered.path)){
      filteredSCE <- importCellRanger(cellRangerDirs = filtered.path, samples = samplename, cellRangerOuts = "/outs/filtered_feature_bc_matrix", class = "Matrix")
    }
} else if(preproc == "SEQC"){
    if(!is.null(droplet.path)){
      dropletSCE <- importSEQC(seqcDirs = droplet.path, samples = samplename, prefix = samplename, class = "Matrix")
    }  
    if(!is.null(filtered.path)){
      filteredSCE <- importSEQC(seqcDirs = filtered.path, samples = samplename, prefix = samplename, class = "Matrix") 
    }
} else if(preproc == "Optimus"){
    if(!is.null(droplet.path)){
        dropletSCE <- importOptimus(OptimusDirs = droplet.path, samples = samplename)
    }
    if(!is.null(filtered.path)){
        filteredSCE <- importOptimus(OptimusDirs = droplet.path, samples = samplename)
    }
} else {
  stop(paste0("'", preproc, "' not supported."))
}

## Read in gene sets for QC
geneSetCollection <- NULL
if(!is.null(gmt)) {
  geneSetCollection <- GSEABase::getGmt(gmt, sep=sep)
}

## Run QC functions
if(!is.null(droplet.path)) {
  message(paste0(date(), " .. Running droplet QC"))    
  dropletSCE <- runDropletQC(sce = dropletSCE)
  
  if(is.null(filtered.path)) {
    ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
    filteredSCE <- dropletSCE[,ix]
  }  
}

message(paste0(date(), " .. Running cell QC"))    
filteredSCE <- runCellQC(sce = filteredSCE, geneSetCollection = geneSetCollection)


## Merge singleCellExperiment objects
if(!is.null(filtered.path) & !is.null(droplet.path)) {
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

if(!is.null(droplet.path)){
  fn <- file.path(directory, samplename, "R", paste0(samplename , "_Droplets.rds"))
  saveRDS(object = mergedDropletSCE, file = fn)
}
fn <- file.path(directory, samplename, "R", paste0(samplename , "_FilteredCells.rds"))
saveRDS(object = mergedFilteredSCE, file = fn)


## ToDo ##
## Export to Python
## Export to flatfile


sessionInfo()
