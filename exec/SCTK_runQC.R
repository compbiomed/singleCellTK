#!/usr/bin/env Rscript --vanilla

##Check to see if necessary packages are installed
#CRAN packages
cran.packages <- c("optparse", "yaml")

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

## Function to parse arguments from yaml file
.parseConfig <- function(sctkConfig, arguments) {
  for (i in seq_along(arguments)) {
    arg <- arguments[i]
    assign(arg, sctkConfig[[arg]], envir = parent.frame())
  }
}

## Check whether python module is available
if (!reticulate::py_module_available(module = "scrublet")) {
    stop("Cannot find python module 'scrublet'. ",
            "Scrublet can be installed on the local machine",
            "with pip (e.g. pip install --user scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
}

##Read in flags from command line using optparse
option_list <- list(optparse::make_option(c("-b", "--base_path"),
        type="character",
        default=NULL,
        help="Base path for the output from the preprocessing algorithm"),
    optparse::make_option(c("-P", "--preproc"),
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
    optparse::make_option(c("-C","--cell_data_path"),
        type="character",
        default=NULL,
        help="The directory contains cell matrix, gene and cell barcodes information. Default is NULL. If 'base_path' is NULL, both 'cell_data_path' and 'raw_data_path' should also be specified."),
    optparse::make_option(c("-R","--raw_data_path"),
        type="character",
        default=NULL,
        help="The directory contains droplet matrix, gene and cell barcodes information. Default is NULL. If 'base_path' is NULL, both 'cell_data_path' and 'raw_data_path' should also be specified."),
    optparse::make_option(c("-S","--split_sample"),
        type="logical",
        default=TRUE,
        help="Save SingleCellExperiment object for each sample. Default is FALSE. If TRUE, all samples will be combined and only one combimed SingleCellExperiment object will be saved."),
    optparse::make_option(c("-r","--raw_data"),
        type="character",
        default=NULL,
        help="The full path of the RDS file or Matrix file of the raw gene count matrix. This would be provided only when --preproc is SceRDS or CountMatrix."),
    optparse::make_option(c("-c","--cell_data"),
        type="character",
        default=NULL,
        help="The full path of the RDS file or Matrix file of the cell count matrix. This would be use only when --preproc is SceRDS or CountMatrix."),
    optparse::make_option(c("-F", "--outputFormat"),
        type="character",
        default=NULL,
        help="The output format of this QC pipeline. Currently, it supports RDS, Flatfile, Python AnnData and HTAN."),
    optparse::make_option(c("-y", "--yamlFile"),
        type="character",
        default=NULL,
        help="YAML file containing parameters called by singleCellTK QC functions. Please check documentation for details."))

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
FilterDir <- opt$cell_data_path 
RawDir <- opt$raw_data_path
Reference <- opt$genome
RawFile <- opt$raw_data
FilterFile <- opt$cell_data
yamlFile <- opt$yamlFile
formats <- opt$outputFormat

if (!is.null(basepath)) { basepath <- unlist(strsplit(opt$base_path, ",")) } 

if (!is.null(FilterDir)) { FilterDir <- unlist(strsplit(opt$cell_data_path, ",")) } 

if (!is.null(RawDir)) { RawDir <- unlist(strsplit(opt$raw_data_path, ",")) } 

if (!is.null(Reference)) { Reference <- unlist(strsplit(opt$genome, ",")) } 

if (!is.null(RawFile)) { RawFile <- unlist(strsplit(opt$raw_data, ",")) }

if (!is.null(FilterFile)) { FilterFile <- unlist(strsplit(opt$cell_data, ",")) } 

if (!is.null(formats)) { formats <- unlist(strsplit(opt$outputFormat, ",")) } 

## Parse parameters for QC algorithms
if (!is.null(yamlFile)) {
    arguments <- c('Params')
    qcParams <- yaml::read_yaml(yamlFile)
    .parseConfig(qcParams, arguments)
} else {
    Params <- list()
}

### checking output formats
if (!all(formats %in% c("R", "Python", "FlatFile", "HTAN"))) {
    warning("Output format must be 'R', 'Python', 'HTAN' or 'FlatFile'. Format ", 
         paste(formats[!format %in% c("R", "Python", "FlatFile", "HTAN")], sep = ","),
         " is not supported now. ") #             "Only output the supported formats in the provided options. "
}

formats <- formats[formats %in% c("R", "Python", "FlatFile", "HTAN")]
message("The output format is [", 
        paste(formats, collapse = ","), "]. ")

if (length(formats) == 0) {
    warning("None of the provided format is supported now. Therefore, the output ", 
        "will be R, Python, FlatFile and HTAN. ")
    formats <- c("R", "Python", "FlatFile", "HTAN")
}

## Checking argument
if (is.null(RawFile) & is.null(RawFile)) {
    if (is.null(basepath)) {
        if ((is.null(FilterDir) || is.null(RawDir))) {
            warning("Both 'cell_data_path' and 'raw_data_path' need to be specified when 'base_path' is NULL.")
        } else {
            # message("'base_path' is NULL. Data is loaded using directories specified by '--cell_data_path' and '--raw_data_path'.")
            if (length(FilterDir) != length(RawDir)) {
                stop("The length of '--cell_data_path' should be the same as the length of '--raw_data_path'.")
            }
            if (length(FilterDir) != length(sample)) {
                stop("The length of '--cell_data_path' should be the same as the length of '--sample'.")
            }
            if (length(FilterDir) != length(process)) {
                stop('The length of "--cell_data_path" should be the same as ',
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
         stop("The length of '--raw_data' and '--cell_data' should be the same when '--preproc' is SceRDS or CountMatrix.")
    }
    if (length(FilterFile) != length(sample)) {
        stop("The length of '--cell_data' should be the same as the length of '--sample'.")
    }
    if (length(FilterFile) != length(process)) {
        stop('The length of "--cell_data" should be the same as ',
                 'the length of "--preproc"!')
    }
}


## Prepare for QC
dropletSCE_list <- list()
cellSCE_list <- list()
geneSetCollection <- NULL
if (!is.null(gmt)) {
    geneSetCollection <- GSEABase::getGmt(gmt, sep=sep)
}

level3Meta <- list()
level4Meta <- list()

for(i in seq_along(process)) {
    preproc <- process[i]
    samplename <- sample[i]
    path <- basepath[i]
    raw <- RawDir[i]
    fil <- FilterDir[i]
    ref <- Reference[i]
    rawFile <- RawFile[i]
    filFile <- FilterFile[i]
    dropletSCE <- NULL
    cellSCE <- NULL

    if (preproc == "BUStools") {
        dropletSCE <- importBUStools(BUStoolsDir = path, sample = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "STARSolo") {
        dropletSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/raw", class = "Matrix", delayedArray=FALSE)
        cellSCE <- importSTARsolo(STARsoloDir = path, sample = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "CellRangerV3") {
        if (!is.null(path)) {
            dropletSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType="raw", class = "Matrix", delayedArray=FALSE)
            cellSCE <- importCellRangerV3(cellRangerDirs = path, sampleNames = samplename, dataType="filtered", class = "Matrix", delayedArray=FALSE)
        } else {
            dropletSCE <- importCellRangerV3Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            cellSCE <- importCellRangerV3Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
        }
    } else if (preproc == "CellRangerV2") {
        if(is.null(ref)){
            stop("The name of genome reference needs to be specified.")
        }
        if (!is.null(path)) {
            dropletSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="raw")
            cellSCE <- importCellRangerV2(cellRangerDirs = path, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="filtered")
        } else {
            dropletSCE <- importCellRangerV2Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            cellSCE <- importCellRangerV2Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
        }
    } else if (preproc == "SEQC") {
        dropletSCE <- importSEQC(seqcDirs = path, samples = samplename, prefix = samplename, class = "Matrix", delayedArray=FALSE)
    } else if (preproc == "Optimus") {
        dropletSCE <- importOptimus(OptimusDirs = path, samples = samplename, delayedArray = FALSE)
        cellSCE <- dropletSCE[,which(dropletSCE$dropletUtils_emptyDrops_IsCell)]
    } else if (preproc == "DropEst") {
        dropletSCE <- importDropEst(sampleDirs=path, dataType="raw", sampleNames=samplename, delayedArray=FALSE)
        cellSCE <- importDropEst(sampleDirs=path, dataType="filtered", sampleNames=samplename, delayedArray=FALSE)
    } else if (preproc == "SceRDS") {
        dropletSCE <- readRDS(rawFile)
        cellSCE <- readRDS(filFile)
    } else if (preproc == "CountMatrix") {
        dropletMM <- data.table::fread(rawFile)
        dropletSCE <- constructSCE(data = dropletMM, samplename = samplename)
        cellMM <- data.table::fread(filFile)
        cellSCE <- constructSCE(data = cellMM, samplename = samplename)
    } else {
        stop(paste0("'", preproc, "' not supported."))
    }
    
    if (!is.null(dropletSCE)) {
        message(paste0(date(), " .. Running droplet QC"))        
        dropletSCE <- runDropletQC(inSCE = dropletSCE, paramsList=Params)
        
        if (is.null(cellSCE)) {
            ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
            cellSCE <- dropletSCE[,ix]
        }    
    }
    
    if (!is.null(cellSCE)) {
        message(paste0(date(), " .. Running cell QC"))        
        cellSCE <- runCellQC(inSCE = cellSCE, geneSetCollection = geneSetCollection, paramsList=Params)
    }
    
    ## merge colData of dropletSCE and FilteredSCE
    mergedDropletSCE <- NULL
    mergedFilteredSCE <- NULL
    if (!is.null(cellSCE) & !is.null(dropletSCE)) {
        mergedDropletSCE <- mergeSCEColData(dropletSCE, cellSCE)
        mergedFilteredSCE <- mergeSCEColData(cellSCE, dropletSCE)
    } else {
        mergedFilteredSCE <- cellSCE
    }

    if (isTRUE(split)) {
        exportSCE(inSCE = mergedDropletSCE, samplename = samplename, directory = directory, type = "Droplets", format=formats)
        exportSCE(inSCE = mergedFilteredSCE, samplename = samplename, directory = directory, type = "Cells", format=formats)
        
        ## Get parameters of QC functions
        getSceParams(inSCE = mergedFilteredSCE, directory = directory, samplename = samplename, writeYAML = TRUE)

        ## generate meta data
        if ("FlatFile" %in% formats) {
            if ("HTAN" %in% formats) {
                meta <- generateMeta(dropletSCE = dropletSCE, cellSCE = cellSCE, samplename = samplename, 
                                    dir = directory, HTAN=TRUE)
            } else {
                meta <- generateMeta(dropletSCE = dropletSCE, cellSCE = cellSCE, samplename = samplename, 
                                    dir = directory, HTAN=FALSE)  
            }

        level3Meta[[i]] <- meta[[1]]
        level4Meta[[i]] <- meta[[2]]

        } else {
            warning("'FlatFile' is not in output format. Skip exporting the manifest file.")
        }
    }

    dropletSCE_list[[samplename]] <- mergedDropletSCE
    cellSCE_list[[samplename]] <- mergedFilteredSCE
}

if (!isTRUE(split)) {
    dropletSCE <- combineSCE(dropletSCE_list)
    cellSCE <- combineSCE(cellSCE_list)

    if (length(sample) > 1) {
        samplename <- paste(sample, collapse="-")
    }

    exportSCE(inSCE = dropletSCE, samplename = samplename, directory = directory, type = "Droplets", format=formats)
    exportSCE(inSCE = cellSCE, samplename = samplename, directory = directory, type = "Cells", format=formats)

    ## Get parameters of QC functions
    getSceParams(inSCE = cellSCE, directory = directory, samplename = samplename, writeYAML = TRUE)

    ## generate meta data
    if ("FlatFile" %in% formats) {
        if ("HTAN" %in% formats) {
            meta <- generateMeta(dropletSCE = dropletSCE, cellSCE = cellSCE, samplename = samplename, 
                                dir = directory, HTAN=TRUE)
        } else {
            meta <- generateMeta(dropletSCE = dropletSCE, cellSCE = cellSCE, samplename = samplename, 
                                dir = directory, HTAN=FALSE)            
        }

        level3Meta <- list(meta[[1]])
        level4Meta <- list(meta[[2]])

    } else {
        warning("'FlatFile' is not in output format. Skip exporting the manifest file.")
    }
}

if ("FlatFile" %in% formats) {
    HTANLevel3 <- do.call(base::rbind, level3Meta)
    HTANLevel4 <- do.call(base::rbind, level4Meta)
    write.csv(HTANLevel3, file = file.path(directory, "level3Meta.csv"))
    write.csv(HTANLevel4, file = file.path(directory, "level4Meta.csv"))
}
