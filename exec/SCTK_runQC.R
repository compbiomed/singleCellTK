#!/usr/bin/env Rscript --vanilla

##Check to see if necessary packages are installed
#CRAN packages
cran.packages <- c("optparse")

cran.package.check <- lapply(cran.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
    }
})

#Bioconductor packages
bioc.packages <- c("singleCellTK", "celda")

bioc.package.check <- lapply(bioc.packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        BiocManager::install(x)
    }
})


#Check which python version is used
Sys.setenv(RETICULATE_PYTHON = '/usr/bin/python3')
reticulate::py_config()


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
        help="Delimiter used in GMT file"),
    optparse::make_option(c("-r","--ref"),
        type="character",
        default=NULL,
        help="The name of genome reference. This is only required for CellRangerV2 data."), 
    optparse::make_option(c("-F","--filtered_expr_path"),
        type="character",
        default=NULL,
        help="The directory contains filtered gene count matrix, gene and cell barcodes information. Default is NULL. If 'base_path' is NULL, both 'filtered_expr_path' and 'raw_expr_path' should also be specified."),
    optparse::make_option(c("-R","--raw_expr_path"),
        type="character",
        default=NULL,
        help="The directory contains raw gene count matrix, gene and cell barcodes information. Default is NULL. If 'base_path' is NULL, both 'filtered_expr_path' and 'raw_expr_path' should also be specified."))

arguments <- optparse::parse_args(optparse::OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options
process <- unlist(strsplit(opt$preproc, ','))
sample <- unlist(strsplit(opt$sample, ','))
directory <- unlist(strsplit(opt$directory, ','))
reference <- unlist(strsplit(opt$ref, ','))
gmt <- opt$gmt
sep <- opt$delim
if (!is.null(opt$base_path)) {
    basepath <- unlist(strsplit(opt$base_path, ','))
} else { 
    basepath <- opt$base_path
}
if (!is.null(opt$filtered_expr_path)) {
    FilterDir <- unlist(strsplit(opt$filtered_expr_path, ','))
} else {
    FilterDir <- opt$filtered_expr_path
}
if (!is.null(opt$raw_expr_path)) {
    RawDir <- unlist(strsplit(opt$raw_expr_path, ','))
} else {
    RawDir <- opt$raw_expr_path
}

## checking argument
if (is.null(basepath)) {
    if (is.null(FilterDir) || is.null(RawDir)) {
        stop("Both 'filtered_expr_path' and 'raw_expr_path' need to be specified when 'base_path' is NULL.")
    } else {
        message("'base_path' is NULL. Data is loaded using directories specified by '--filtered_expr_path' and '--raw_expr_path'.")
        if (length(FilterDir) != length(RawDir)) {
            stop("The length of '--filtered_expr_path' should be the same as the length of '--raw_expr_path'.")
        }
        if (length(FilterDir) != length(sample)) {
            stop("The length of '--filtered_expr_path' should be the same as the length of '--sample'.")
        }
        if (length(FilterDir) != length(process)) {
            stop('The length of "--filtered_expr_path" should be the same as ',
                     'the length of "--preproc"!')
        }
    }
} else {
    if (length(basepath) != length(process)) {
        stop('The length of "--base_path" should be the same as ',
                 'the length of "--preproc"!')
    }
    if (length(basepath) != length(sample)) {
        stop('The length of "--base_path" should be the same as ',
                 'the length of "--sample"!')    
    }
}

if (length(reference) != sum(process == 'CellRangerV2')) {
    stop('The length of "--ref" should be the same as ',
             'the number of "CellRangerV2" in the "--preproc"!')        
}

dropletSCE_list <- list()
filteredSCE_list <- list()
geneSetCollection <- NULL
if (!is.null(gmt)) {
    geneSetCollection <- GSEABase::getGmt(gmt, sep=sep)
}

for(i in 1:length(process)) {
    preproc <- process[i]
    samplename <- sample[i]
    path <- basepath[i]
    raw <- RawDir[i]
    fil <- FilterDir[i]
    ref <- reference[i]
    dropletSCE <- NULL
    filteredSCE <- NULL
    
    if (preproc == "BUStools") {
        dropletSCE <- importBUStools(BUStoolsDir = path, sample = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "STARSolo") {
        dropletSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/raw", class = "Matrix", delayedArray=FALSE)
        filteredSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix")
    } else if (preproc == "CellRangerV3") {
        if (!is.null(path)) {
            dropletSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType='raw', class = 'Matrix', delayedArray=FALSE)
            filteredSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType='filtered', class = "Matrix")
        } else {
            dropletSCE <- importCellRangerV3Sample(dataDir = raw, sampleName = samplename, class = 'Matrix', delayedArray=FALSE)
            filteredSCE <- importCellRangerV3Sample(dataDir = fil, sampleName = samplename, class = 'Matrix', delayedArray=TRUE)
        }
    } else if (preproc == "CellRangerV2") {
        if(is.null(ref)){
            stop("The name of genome reference needs to be specified.")
        } else {
            rawOuts <- paste0("outs/raw_gene_bc_matrices/", ref)
            filterOuts <- paste0("outs/filtered_gene_bc_matrices/", ref)
        }
        dropletSCE <- importCellRanger(cellRangerDirs = path, sampleNames = samplename, gzipped = FALSE, cellRangerOuts = rawOuts, class = "Matrix", matrixFileNames = "matrix.mtx", featuresFileNames = "genes.tsv", barcodesFileNames = "barcodes.tsv", delayedArray=FALSE)
        filteredSCE <- importCellRanger(cellRangerDirs = path, sampleNames = samplename, gzipped = FALSE, cellRangerOuts = filterOuts, class = "Matrix", matrixFileNames = "matrix.mtx", featuresFileNames = "genes.tsv", barcodesFileNames = "barcodes.tsv", delayedArray=FALSE)

        # if (!is.null(path)) {
        #     dropletSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class='Matrix', delayedArray = FALSE, reference = ref, dataTypeV2='raw')
        #     filteredSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class='Matrix', delayedArray = FALSE, reference = ref, dataTypeV2='filtered')
        # } else {
        #     dropletSCE <- importCellRangerV2Sample(dataDir = raw, sampleName = samplename, class = 'Matrix', delayedArray=FALSE)
        #     filteredSCE <- importCellRangerV2Sample(dataDir = fil, sampleName = samplename, class = 'Matrix', delayedArray=TRUE)
        # }
    } else if (preproc == "SEQC") {
        dropletSCE <- importSEQC(seqcDirs = path, samples = samplename, prefix = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "Optimus") {
        dropletSCE <- importOptimus(OptimusDirs = path, samples = samplename, delayedArray = FALSE)
        filteredSCE <- dropletSCE[,which(dropletSCE$dropletUtils_emptyDrops_IsCell)]
    } else {
        stop(paste0("'", preproc, "' not supported."))
    }
    
    if (!is.null(dropletSCE)) {
        message(paste0(date(), " .. Running droplet QC"))        
        dropletSCE <- runDropletQC(inSCE = dropletSCE) # "emptyDrops", 
        
        if (is.null(filteredSCE)) {
            ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
            filteredSCE <- dropletSCE[,ix]
        }    
    }
    
    if (!is.null(filteredSCE)) {
        message(paste0(date(), " .. Running cell QC"))        
        filteredSCE <- runCellQC(inSCE = filteredSCE, geneSetCollection = geneSetCollection)
    }
    
    dropletSCE_list[[samplename]] <- dropletSCE
    filteredSCE_list[[samplename]] <- filteredSCE
}

dropletSCE <- do.call(BiocGenerics::cbind, dropletSCE_list)
filteredSCE <- do.call(BiocGenerics::cbind, filteredSCE_list)

mergedDropletSCE <- NULL
mergedFilteredSCE <- NULL
if (!is.null(filteredSCE) & !is.null(dropletSCE)) {
    mergedDropletSCE <- mergeSCEColData(dropletSCE, filteredSCE)
    mergedFilteredSCE <- mergeSCEColData(filteredSCE, dropletSCE)
} else {
    mergedFilteredSCE <- filteredSCE
}

if (length(sample) > 1) {
    samplename <- paste(sample, collapse='_')
}
## Create directories and save objects
dir.create(file.path(directory, samplename), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "R"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "Python"), showWarnings = TRUE, recursive = TRUE)
dir.create(file.path(directory, samplename, "FlatFile"), showWarnings = TRUE, recursive = TRUE)

if (!is.null(mergedDropletSCE)) {
    ## Export to R 
    fn <- file.path(directory, samplename, "R", paste0(samplename , "_Droplets.rds"))
    saveRDS(object = mergedDropletSCE, file = fn)
    
    ## Export to flatfile
    fn <- file.path(directory, samplename, "FlatFile", "Droplets")
    exportSCEtoFlatFile(mergedDropletSCE, outputDir = fn)
}
if (!is.null(mergedFilteredSCE)) {
    ## Export to R        
    fn <- file.path(directory, samplename, "R", paste0(samplename , "_FilteredCells.rds"))
    saveRDS(object = mergedFilteredSCE, file = fn)
    
    ## Export to flatfile    
    fn <- file.path(directory, samplename, "FlatFile", "FilteredCells")
    exportSCEtoFlatFile(mergedFilteredSCE, outputDir = fn)
}    