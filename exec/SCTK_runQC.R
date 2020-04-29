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

#Output function 
sceOutput <- function(dropletSCE, filteredSCE, samplename, directory){
    mergedDropletSCE <- NULL
    mergedFilteredSCE <- NULL
    if (!is.null(filteredSCE) & !is.null(dropletSCE)) {
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
  
    if (!is.null(mergedDropletSCE)) {
        ## Export to R 
        fn <- file.path(directory, samplename, "R", paste0(samplename , "_Droplets.rds"))
        saveRDS(object = mergedDropletSCE, file = fn)
    
        ## Export to flatfile
        fn <- file.path(directory, samplename, "FlatFile", "Droplets")
        exportSCEtoFlatFile(mergedDropletSCE, outputDir = fn)
    
        ## Export to Python AnnData
        fn <- file.path(directory, samplename, "Python", "Droplets")
        exportSCEtoAnnData(mergedDropletSCE, outputDir=fn, compression='gzip', sample=samplename)
    }

    if (!is.null(mergedFilteredSCE)) {
        ## Export to R        
        fn <- file.path(directory, samplename, "R", paste0(samplename , "_FilteredCells.rds"))
        saveRDS(object = mergedFilteredSCE, file = fn)
        
        ## Export to flatfile    
        fn <- file.path(directory, samplename, "FlatFile", "FilteredCells")
        exportSCEtoFlatFile(mergedFilteredSCE, outputDir = fn)
        
        ## Export to Python AnnData
        fn <- file.path(directory, samplename, "Python", "FilteredCells")
        exportSCEtoAnnData(mergedFilteredSCE, outputDir=fn, compression='gzip', sample=samplename)
  }
}

# Function of combining SingleCellExperiment object
combineSCE <- function(sceList){
    qcList <- sapply(sceList, function(x) {colnames(x@colData)})
    qcMetNum <- sapply(qcList, length)
  
    if (var(qcMetNum) != 1) { ##some QC alrorithms failed for some samples
        qcMetrics <- base::Reduce(union, qcList)
    
    for (i in seq_along(sceList)) {
        sce <- sceList[[i]]
        missQC <- qcMetrics[!qcMetrics %in% colnames(sce@colData)]
  
        if (length(missQC) != 0) {
            missColDat <- S4Vectors::DataFrame(sapply(missQC, function(x){rep(NA, ncol(sce))}))
            colData(sce) <- cbind(colData(sce), missColDat)
            sceList[[i]] <- sce      
        }
    }
  }
    sce <- do.call(BiocGenerics::cbind, sceList)
    return(sce)
}

## create SCE from csv or txt input
constructSCE <- function(data, samplename){
    gene <- data[[1]]
    data <- data[, -1]
    barcode <- colnames(data)
    mat <- methods::as(data, "Matrix")
    dimnames(mat) <- list(gene, barcode)
    coln <- paste(samplename, barcode, sep = '_')
  
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = mat))
    SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(feature = gene)
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(barcode,
        column_name = coln,
        sample = samplename,
        row.names = coln)
  
    return(sce)
}

##Read in flags from command line using optparse
option_list <- list(optparse::make_option(c("-b", "--base_path"),
        type="character",
        default=NULL,
        help="Base path for the output from the preprocessing algorithm"),
    optparse::make_option(c("-p", "--preproc"),
        type = "character",
        default="CellRangerV3",
        help="Algorithm used for preprocessing. One of 'CellRangerV2', 'CellRangerV3', 'BUStools', 'STARSolo', 'SEQC', 'Optimus', 'DropEst', 'SceRDS', 'CountMatrix'"),
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
    optparse::make_option(c("-G","--genome"),
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
        help="The directory contains raw gene count matrix, gene and cell barcodes information. Default is NULL. If 'base_path' is NULL, both 'filtered_expr_path' and 'raw_expr_path' should also be specified."),
    optparse::make_option(c("-S","--split_sample"),
        type="logical",
        default=FALSE,
        help="Save SingleCellExperiment object for each sample. Default is FALSE. If TRUE, all samples will be combined and only one combimed SingleCellExperiment object will be saved."),
    optparse::make_option(c("-r","--raw_data"),
        type="character",
        default=NULL,
        help="The full path of the RDS file or Matrix file of the raw gene count matrix. This would be provided only when --preproc is SceRDS or CountMatrix."),
    optparse::make_option(c("-f","--filtered_data"),
        type="character",
        default=NULL,
        help="The full path of the RDS file or Matrix file of the filtered gene count matrix. This would be use only when --preproc is SceRDS or CountMatrix."))

## Define arguments
arguments <- optparse::parse_args(optparse::OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options
process <- unlist(strsplit(opt$preproc, ","))
sample <- unlist(strsplit(opt$sample, ","))
directory <- unlist(strsplit(opt$directory, ","))
gmt <- opt$gmt
sep <- opt$delim
split <- opt$split_sample
basepath <- opt$base_path
FilterDir <- opt$filtered_expr_path 
RawDir <- opt$raw_expr_path
Reference <- opt$genome
RawFile <- opt$raw_data
FilterFile <- opt$filtered_data

if (!is.null(basepath)) { basepath <- unlist(strsplit(opt$base_path, ",")) } 

if (!is.null(FilterDir)) { FilterDir <- unlist(strsplit(opt$filtered_expr_path, ",")) } 

if (!is.null(RawDir)) { RawDir <- unlist(strsplit(opt$raw_expr_path, ",")) } 

if (!is.null(Reference)) { Reference <- unlist(strsplit(opt$genome, ",")) } 

if (!is.null(RawFile)) { RawFile <- unlist(strsplit(opt$raw_data, ",")) }

if (!is.null(FilterFile)) { FilterFile <- unlist(strsplit(opt$filtered_data, ",")) } 

## checking argument
if (is.null(RawFile) & is.null(RawFile)) {
    if (is.null(basepath)) {
        if ((is.null(FilterDir) || is.null(RawDir))) {
            warning("Both 'filtered_expr_path' and 'raw_expr_path' need to be specified when 'base_path' is NULL.")
        } else {
            # message("'base_path' is NULL. Data is loaded using directories specified by '--filtered_expr_path' and '--raw_expr_path'.")
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

    if (length(Reference) != sum(process == 'CellRangerV2')) {
        stop('The length of "--ref" should be the same as ',
                 'the number of "CellRangerV2" in the "--preproc"!')        
    }        
}

if (!is.null(RawFile) | !is.null(FilterFile)) {
    if (length(RawFile) != length(FilterFile)) {
         stop("The length of '--raw_data' and '--filtered_data' should be the same when '--preproc' is SceRDS or CountMatrix.")
    }
    if (length(FilterFile) != length(sample)) {
        stop("The length of '--filtered_data' should be the same as the length of '--sample'.")
    }
    if (length(FilterFile) != length(process)) {
        stop('The length of "--filtered_data" should be the same as ',
                 'the length of "--preproc"!')
    }
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
    ref <- Reference[i]
    rawFile <- RawFile[i]
    filFile <- FilterFile[i]
    dropletSCE <- NULL
    filteredSCE <- NULL

    if (preproc == "BUStools") {
        dropletSCE <- importBUStools(BUStoolsDir = path, sample = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "STARSolo") {
        dropletSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/raw", class = "Matrix", delayedArray=FALSE)
        filteredSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "CellRangerV3") {
        if (!is.null(path)) {
            dropletSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType="raw", class = "Matrix", delayedArray=FALSE)
            filteredSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType="filtered", class = "Matrix", delayedArray=FALSE)
        } else {
            dropletSCE <- importCellRangerV3Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            filteredSCE <- importCellRangerV3Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
        }
    } else if (preproc == "CellRangerV2") {
        if(is.null(ref)){
            stop("The name of genome reference needs to be specified.")
        }
        if (!is.null(path)) {
            dropletSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="raw")
            filteredSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="filtered")
        } else {
            dropletSCE <- importCellRangerV2Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            filteredSCE <- importCellRangerV2Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
        }
    } else if (preproc == "SEQC") {
        dropletSCE <- importSEQC(seqcDirs = path, samples = samplename, prefix = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "Optimus") {
        dropletSCE <- importOptimus(OptimusDirs = path, samples = samplename, delayedArray = FALSE)
        filteredSCE <- dropletSCE[,which(dropletSCE$dropletUtils_emptyDrops_IsCell)]
    } else if (preproc == "DropEst") {
        dropletSCE <- importDropEst(sampleDirs=path, dataType="raw", sampleNames=samplename, delayedArray=FALSE)
        filteredSCE <- importDropEst(sampleDirs=path, dataType="filtered", sampleNames=samplename, delayedArray=FALSE)
    } else if (preproc == "SceRDS") {
        dropletSCE <- readRDS(rawFile)
        filteredSCE <- readRDS(filFile)
    } else if (preproc == "CountMatrix") {
        dropletMM <- data.table::fread(rawFile)
        dropletSCE <- constructSCE(data = dropletMM, samplename = samplename)
        filteredMM <- data.table::fread(filFile)
        filteredSCE <- constructSCE(data = filteredMM, samplename = samplename)
    } else {
        stop(paste0("'", preproc, "' not supported."))
    }
    
    if (!is.null(dropletSCE)) {
        message(paste0(date(), " .. Running droplet QC"))        
        dropletSCE <- runDropletQC(inSCE = dropletSCE)
        
        if (is.null(filteredSCE)) {
            ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
            filteredSCE <- dropletSCE[,ix]
        }    
    }
    
    if (!is.null(filteredSCE)) {
        message(paste0(date(), " .. Running cell QC"))        
        filteredSCE <- runCellQC(inSCE = filteredSCE, geneSetCollection = geneSetCollection)
    }
    
    if (isTRUE(split)) {
        sceOutput(dropletSCE=dropletSCE, filteredSCE=filteredSCE, samplename=samplename, directory=directory)
    }

    dropletSCE_list[[samplename]] <- dropletSCE
    filteredSCE_list[[samplename]] <- filteredSCE
}

if (!isTRUE(split)){
    dropletSCE <- combineSCE(dropletSCE_list)
    filteredSCE <- combineSCE(filteredSCE_list)

    if (length(sample) > 1) {
        samplename <- "Combined"
    }

    sceOutput(dropletSCE=dropletSCE, filteredSCE=filteredSCE, samplename=samplename, directory=directory)
}
