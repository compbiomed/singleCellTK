#' Export data in SingleCellExperiment object
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object 
#' that contains the data. QC metrics are stored in colData of the 
#' singleCellExperiment object.
#' @param samplename Sample name. This will be used as name of subdirectories
#' and the prefix of flat file output. Default is 'sample'.
#' @param type Type of data. The type of data stored in SingleCellExperiment object. 
#' It can be 'Droplets'(raw droplets matrix) or 'Cells' (cells matrix). 
#' @param directory Output directory. Default is './'.
#' @param format The format of output. It currently supports flat files, rds files
#' and python h5 files. It can output multiple formats. Default: c("SCE", "AnnData", "FlatFile", "HTAN"). 
#' @export
exportSCE <- function(inSCE, 
                      samplename = "sample", 
                      directory = "./", 
                      type = "Cells",
                      format = c("SCE", "AnnData", "FlatFile", "HTAN")) {
  
    if (any(!format %in% c("SCE", "AnnData", "FlatFile", "HTAN"))) {
        warning("Output format must be 'SCE', 'AnnData', 'HTAN' or 'FlatFile'. Format ", 
             paste(format[!format %in% c("SCE", "AnnData", "FlatFile", "HTAN")], sep = ","),
             " is not supported now. ") #             "Only output the supported formats in the provided options. "
    }

    format <- format[format %in% c("SCE", "AnnData", "FlatFile", "HTAN")]
    message("The output format is [", 
            paste(format, collapse = ","), "]. ")

    if (length(format) == 0) {
        warning("None of the provided format is supported now. Therefore, the output ", 
            "will be SCE, AnnData, FlatFile and HTAN. ")
        format <- c("SCE", "AnnData", "FlatFile", "HTAN")
    }

    ## Create directories and save objects
    dir.create(file.path(directory, samplename), showWarnings = TRUE, recursive = TRUE)
  
    if ("SCE" %in% format) {
        ## Export to R
        fp <- file.path(directory, samplename, "R")
        dir.create(fp, showWarnings = TRUE, recursive = TRUE)
        fn <- file.path(fp, paste0(samplename , paste0("_", type, ".rds"))) 
        saveRDS(object = inSCE, file = fn)
    } 

    if ("FlatFile" %in% format) {
        ## Export to flatfile
        fp <- file.path(directory, samplename, "FlatFile")
        dir.create(fp, showWarnings = TRUE, recursive = TRUE)
        fn <- file.path(fp, type)
        exportSCEtoFlatFile(inSCE, outputDir = fn, sample=samplename)
    }

    if ("AnnData" %in% format) {
        ## Export to Python AnnData
        fp <- file.path(directory, samplename, "Python")
        dir.create(fp, showWarnings = TRUE, recursive = TRUE)
        fn <- file.path(fp, type)
        exportSCEtoAnnData(inSCE, outputDir=fn, compression='gzip', prefix=samplename)
    }
}


#' Combine a list of SingleCellExperiment objects as one SingleCellExperiment object
#' @param sceList A list contains \link[SingleCellExperiment]{SingleCellExperiment} objects 
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object which combines all 
#' objects in sceList. The colData is merged.  
#' @export
#' @importFrom SummarizedExperiment colData colData<-
combineSCE <- function(sceList) {
    qcList <- sapply(sceList, function(x) {colnames(x@colData)})
    qcMetNum <- sapply(qcList, length)
  
    if (stats::var(qcMetNum) != 1) { ##some QC alrorithms failed for some samples
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

#' Generate manifest file for droplet and cell count data
#' @param dropletSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' droplet count matrix data
#' @param cellSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' cell count matrix data
#' @param samplename The sample name of the \link[SingleCellExperiment]{SingleCellExperiment} objects
#' @param dir The output directory of the SCTK QC pipeline. 
#' @param HTAN Whether generate manifest file with the format required by HTAN. Default is TRUE. 
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object which combines all 
#' objects in sceList. The colData is merged.  
#' @export
#' @importFrom SummarizedExperiment assay colData
generateMeta <- function(dropletSCE, 
                          cellSCE, 
                          samplename, 
                          dir,
                          HTAN=TRUE) {
  level3List <- list()
  level4List <- list()
  
  directory <- file.path(basename(dir), samplename)
  filterDir <- file.path(directory, 'FlatFile', 'Cells')
  rawDir <- file.path(directory, 'FlatFile', 'Droplets')
  
  
  WorkFlowData = c(
    WorkFlow = 'singleCellTK QC pipeline', 
    WorkFlowVer = paste('singleCellTK', utils::sessionInfo()$otherPkgs$singleCellTK$Version, sep=':'),
    ParRaw = 'Ran perCellQC, EmptyDrops and barcodeRankDrops using singleCellTK',
    ParFiltered = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    PardecontX = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    ColData = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    DecontXUMAP = 'UMAP dimension reduction generated by decontX',
    ScrubletTSNE = 'tSNE dimension reduction generated by Scrublet',
    ScrubletUMAP = 'UMAP dimension reduction generated by Scrublet'
  )
  
  data <- list(
    'Raw' = c(CellNum = ncol(dropletSCE),
              MedianReads = stats::median(colData(dropletSCE)$sum),
              MedianGenes = stats::median(colData(dropletSCE)$detected),
              DataType = 'Droplet Matrix',
              FileName = file.path(rawDir, 'assays', paste0(samplename,'_counts.mtx.gz'))), 
    
    'decontX' = c(CellNum = ncol(cellSCE),
                  MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts'))),
                  MedianGenes = stats::median(apply(assay(cellSCE, 'decontXcounts'), 2, function(x){sum(x>0)})),
                  DataType = 'Decontaminated cell matrix return returned by runDecontX',
                  FileName = file.path(filterDir, 'assays', paste0(samplename,'_decontXcounts.mtx.gz'))), 
    
    'Filtered' = c(CellNum = ncol(cellSCE),
                   MedianReads = stats::median(colData(cellSCE)$sum),
                   MedianGenes = stats::median(colData(cellSCE)$detected),
                   DataType = 'Cell Matrix', 
                   FileName = file.path(filterDir, 'assays', paste0(samplename,'_counts.mtx.gz')),
                   ColData = file.path(filterDir, paste0(samplename,'_colData.txt.gz')),
                   DecontXUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP.txt.gz')),
                   ScrubletTSNE = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_TSNE.txt.gz')),
                   ScrubletUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_TSNE.txt.gz'))
    ))

  for (type in c('Raw', 'Filtered', 'decontX')) {
    level3List[[type]] <- data.frame(
      SAMPLE = samplename, DATA_TYPE = data[[type]]['DataType'],
      MATRIX_TYPE = 'Raw Counts', DATA_CATEGORY = 'Gene Expression',
      CELL_MEDIAN_NUM_READS = data[[type]]['MedianReads'],
      CELL_MEDIAN_NUM_GENES = data[[type]]['MedianGenes'],
      CELL_TOTAL = data[[type]]['CellNum'],
      FILE_NAME = data[[type]]['FileName'],
      WORKFLOW_TYPE = WorkFlowData['WorkFlow'],
      WORKFLOW_PARAMETERS = WorkFlowData[paste0('Par', type)],
      WORKFLOW_VERSION = WorkFlowData['WorkFlowVer'],
      stringsAsFactors = FALSE)
    
    if (isTRUE(HTAN)) {
      level3List[[type]] <- cbind(level3List[[type]],HTAN_BIOSPECIMEN_ID = samplename,
                                  HTAN_PARENT_ID = '', HTAN_PARENT_FILE_ID = '')
    } 
       
    if (type == 'Filtered') {
      for (metric in c('ColData', 'DecontXUMAP', 'ScrubletTSNE', 'ScrubletUMAP')) {
        level4List[[metric]] <- data.frame(   
          SAMPLE = samplename, FILE_NAME = data[[type]][metric], 
          WORKFLOW_TYPE = WorkFlowData[metric],
          WORKFLOW_PARAMETERS = file.path(directory, paste0(samplename, '_QCParameters.yaml')),
          WORKFLOW_VERSION = WorkFlowData['WorkFlowVer'])
        if (isTRUE(HTAN)) {
          level4List[[metric]] <- cbind(level4List[[metric]],  HTAN_BIOSPECIMEN_ID = samplename,
                                        HTAN_PARENT_ID = '',HTAN_PARENT_FILE_ID = data[[type]]['FileName'])
        }
      }
    }
  }
  
  level3Meta <- do.call(base::rbind, level3List)
  level4Meta <- do.call(base::rbind, level4List)
  return(list(level3Meta, level4Meta))
}

#' Extract QC parameters from the SingleCellExperiment object
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object. 
#' @param skip Skip extracting the parameters of the provided QC functions. 
#' @param ignore Skip extracting the content within QC functions.
#' @param directory The output directory of the SCTK_runQC.R pipeline.
#' @param samplename The sample name of the \link[SingleCellExperiment]{SingleCellExperiment} objects.
#' @param writeYAML Whether output yaml file to store parameters. Default if TRUE. If FALSE, 
#' return character object.   
#' @export
getSceParams <- function(inSCE, 
                         skip = c("scrublet", "runDecontX"), 
                         ignore = c("algorithms", "estimates","contamination",
                                    "z","sample","rank","BPPARAM","batch","geneSetCollection"), 
                         directory = './', 
                         samplename = '',
                         writeYAML = TRUE) {
  
  meta <- S4Vectors::metadata(inSCE)
  algos <- names(meta)[!names(meta) %in% skip]
  outputs <- '---'
  parList <- list()
  dir <- file.path(directory, samplename)
  
  for (algo in algos) {
    params <- meta[[algo]]
    if (length(params) == 1) {params <- params[[1]]} ### extract params from sublist
    params <- params[which(!names(params) %in% ignore)]
    parList[[algo]] <- params
  }

  outputs <- paste(outputs, yaml::as.yaml(parList), sep='\n')
  if (isTRUE(writeYAML)) {
    filename <- paste0(samplename, '_QCParameters.yaml')
    cat(outputs, file=file.path(dir, filename))
  } else {
    return(outputs)
  }
}


#' Create SingleCellExperiment object from csv or txt input
#' @param data A \link[data.table]{data.table} object containing the count matrix. 
#' @param samplename The sample name of the data.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the count matrix.  
#' @export
constructSCE <- function(data, samplename) {
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

#' Create SingleCellExperiment object from command line input arguments
#' @param preproc Method used to preprocess the data. It's one of the path provided in --preproc argument.
#' @param path Base path of the dataset. It's one of the path provided in --bash_path argument.  
#' @param samplename The sample name of the data. It's one of the path provided in --sample argument.  
#' @param raw The directory contains droplet matrix, gene and cell barcodes information. It's one of the path provided in --raw_data_path argument.
#' @param fil The directory contains cell matrix, gene and cell barcodes information. It's one of the path provided in --cell_data_path argument.
#' @param ref The name of reference used by cellranger. Only need for CellrangerV2 data. 
#' @param rawFile The full path of the RDS file or Matrix file of the raw gene count matrix. It's one of the path provided in --raw_data argument.
#' @param filFile The full path of the RDS file or Matrix file of the cell count matrix. It's one of the path provided in --cell_data argument.
#' @param dataType Type of the input. It can be "Both", "Droplet" or "Cell". It's one of the path provided in --genome argument.
#' @return A list of \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the droplet or cell data or both,depending on the dataType that users provided.  
#' @export
qcInputProcess <- function(preproc,
    samplename,
    path,
    raw,
    fil,
    ref,
    rawFile,
    filFile,
    dataType) {

    dropletSCE <- NULL
    cellSCE <- NULL

    if (preproc == "BUStools") {
        dropletSCE <- importBUStools(BUStoolsDirs = path, samples = samplename, class = "Matrix", delayedArray=FALSE)
        return(list(dropletSCE, cellSCE))
    } 

    if (preproc == "SEQC") {
        dropletSCE <- importSEQC(seqcDirs = path, samples = samplename, prefix = samplename, class = "Matrix", delayedArray=FALSE)
        return(list(dropletSCE, cellSCE))
    }

    if (preproc == "STARSolo") {
        if (dataType == "Both") {
            dropletSCE <- importSTARsolo(STARsoloDirs = path, samples = samplename, STARsoloOuts = "Gene/raw", class = "Matrix", delayedArray=FALSE)
            cellSCE <- importSTARsolo(STARsoloDirs = path, samples = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix", delayedArray=FALSE)
        } else if (dataType == "Cell") {
            cellSCE <- importSTARsolo(STARsoloDirs = path, samples = samplename, STARsoloOuts = "Gene/filtered", class = "Matrix", delayedArray=FALSE)
        } else if (dataType == "Droplet") {
            dropletSCE <- importSTARsolo(STARsoloDirs = path, samples = samplename, STARsoloOuts = "Gene/raw", class = "Matrix", delayedArray=FALSE)
        }
        return(list(dropletSCE, cellSCE))
    } 

    if (preproc == "CellRangerV3") {
        if (!is.null(path)) {
            if (dataType == "Both") {
                dropletSCE <- importCellRangerV3(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, dataType="raw", class = "Matrix", delayedArray=FALSE)
                cellSCE <- importCellRangerV3(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, dataType="filtered", class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Cell") {
                cellSCE <- importCellRangerV3(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, dataType="filtered", class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Droplet") {
                dropletSCE <- importCellRangerV3(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, dataType="raw", class = "Matrix", delayedArray=FALSE)
            }
        } else {
            if (dataType == "Both") {
                dropletSCE <- importCellRangerV3Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
                cellSCE <- importCellRangerV3Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Cell") {
                cellSCE <- importCellRangerV3Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Droplet") {
                dropletSCE <- importCellRangerV3Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            }
        }
        return(list(dropletSCE, cellSCE))
    }

    if (preproc == "CellRangerV2") {
        if (!is.null(path)) {
            if (dataType == "Both") {
                dropletSCE <- importCellRangerV2(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="raw")
                cellSCE <- importCellRangerV2(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="filtered")
            } else if (dataType == "Cell") {
                cellSCE <- importCellRangerV2(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="filtered")
            } else if (dataType == "Droplet") {
                dropletSCE <- importCellRangerV2(cellRangerDirs = path, sampleDirs = samplename, sampleNames = samplename, class="Matrix", delayedArray = FALSE, reference = ref, dataTypeV2="raw")
            }
        } else {
            if (dataType == "Both") {
                dropletSCE <- importCellRangerV2Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
                cellSCE <- importCellRangerV2Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Cell") {
                cellSCE <- importCellRangerV2Sample(dataDir = fil, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            } else if (dataType == "Droplet") {
                dropletSCE <- importCellRangerV2Sample(dataDir = raw, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
            }
        }
        return(list(dropletSCE, cellSCE))
    }

    if (preproc == "Optimus") {
        ## by default has both droplet and cell data. 
        dropletSCE <- importOptimus(OptimusDirs = path, samples = samplename, delayedArray = FALSE)
        cellSCE <- dropletSCE[,which(dropletSCE$dropletUtils_emptyDrops_IsCell)]

        if (dataType == "Cell") {
            dropletSCE <- NULL
        } else if (dataType == "Droplet") {
            cellSCE <- NULL
        }
        return(list(dropletSCE, cellSCE))
    }

    if (preproc == "DropEst") {
        if (dataType == "Both") {
            dropletSCE <- importDropEst(sampleDirs=path, dataType="raw", sampleNames=samplename, delayedArray=FALSE)
            cellSCE <- importDropEst(sampleDirs=path, dataType="filtered", sampleNames=samplename, delayedArray=FALSE)
        } else if (dataType == "Cell") {
            cellSCE <- importDropEst(sampleDirs=path, dataType="filtered", sampleNames=samplename, delayedArray=FALSE)
        } else if (dataType == "Droplet") {
            dropletSCE <- importDropEst(sampleDirs=path, dataType="raw", sampleNames=samplename, delayedArray=FALSE)
        }
        return(list(dropletSCE, cellSCE))    
    }

    if (preproc == "SceRDS") {
        if (dataType == "Both") {
            dropletSCE <- readRDS(rawFile)
            cellSCE <- readRDS(filFile)
        } else if (dataType == "Cell") {
            cellSCE <- readRDS(filFile)
        } else if (dataType == "Droplet") {
            dropletSCE <- readRDS(rawFile)
        }
        return(list(dropletSCE, cellSCE))
    }

    if (preproc == "CountMatrix") {
        if (dataType == "Both") {
            dropletMM <- data.table::fread(rawFile)
            dropletSCE <- constructSCE(data = dropletMM, samplename = samplename)
            cellMM <- data.table::fread(filFile)
            cellSCE <- constructSCE(data = cellMM, samplename = samplename)
        } else if (dataType == "Cell") {
            cellMM <- data.table::fread(filFile)
            cellSCE <- constructSCE(data = cellMM, samplename = samplename)
        } else if (dataType == "Droplet") {
            dropletMM <- data.table::fread(rawFile)
            dropletSCE <- constructSCE(data = dropletMM, samplename = samplename)
        }
        return(list(dropletSCE, cellSCE))
    }

    ## preproc is not one of the method above. Stop the pipeline. 
    stop(paste0("'", preproc, "' not supported."))
}