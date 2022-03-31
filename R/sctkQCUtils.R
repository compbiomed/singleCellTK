.convert_hex_to_int <- function(hex) {
  code = paste0("hash = str(int('", hex, "', base=16))[-10:]")
  reticulate::py_run_string(code)
  reticulate::py$hash
}

#' Export data in Seurat object
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object
#' that contains the data. QC metrics are stored in colData of the
#' singleCellExperiment object.
#' @param outputDir Path to the directory where outputs will be written. Default is the current working directory.
#' @param prefix Prefix to use for the name of the output file. Default \code{"sample"}.
#' @param overwrite Boolean. Whether overwrite the output if it already exists in the outputDir. Default \code{TRUE}.
#' @param copyColData Boolean. Whether copy 'colData' of SCE object to the 'meta.data' of Seurat object. Default \code{TRUE}.
#' @param copyReducedDim Boolean. Whether copy 'reducedDims' of the SCE object to the 'reductions' of Seurat object. Default \code{TRUE}.
#' @param copyDecontX Boolean. Whether copy 'decontXcounts' assay of the SCE object to the 'assays' of Seurat object. Default \code{TRUE}.
#' @return Generates a Seurat object containing data from \code{inSCE}.
#' @export

exportSCEToSeurat <- function(inSCE, prefix="sample", outputDir="./", overwrite=TRUE,
                              copyColData=TRUE, copyReducedDim=TRUE,
                              copyDecontX=TRUE) {
  Seurat <- singleCellTK::convertSCEToSeurat(inSCE, countsAssay = "counts", 
                                             copyColData = copyColData,
                                             copyReducedDim = copyReducedDim,
                                             copyDecontX = copyDecontX)

  fileName <- paste0(prefix,"_Seurat.RDS")
  filePath <- file.path(outputDir,fileName)

  if (!dir.exists(outputDir)) {
    message("outputDir does not exists. Create the directory. ")
    dir.create(outputDir)
  }

  if (file.exists(filePath) && !isTRUE(overwrite)) {
    stop(paste0(path, " already exists. Change 'outputDir' or set 'overwrite' to TRUE."))
  }


  saveRDS(Seurat, filePath)
}

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
#' @return Generates a file containing data from \code{inSCE}, in specified \code{format}.
#' @examples
#' data(scExample)
#' \dontrun{
#' exportSCE(sce, format = "SCE")
#' }
#' @export
exportSCE <- function(inSCE,
                      samplename = "sample",
                      directory = "./",
                      type = "Cells",
                      format = c("SCE", "AnnData", "FlatFile", "HTAN", "Seurat")) {

    if (any(!format %in% c("SCE", "AnnData", "FlatFile", "HTAN", "Seurat"))) {
        warning("Output format must be 'SCE', 'AnnData', 'HTAN', 'Seurat' or 'FlatFile'. Format ",
             paste(format[!format %in% c("SCE", "AnnData", "FlatFile", "HTAN")], sep = ","),
             " is not supported now. ") #             "Only output the supported formats in the provided options. "
    }

    format <- format[format %in% c("SCE", "AnnData", "FlatFile", "HTAN", "Seurat")]
    message("The output format is [",
            paste(format, collapse = ","), "]. ")

    if (length(format) == 0) {
        warning("None of the provided format is supported now. Therefore, the output ",
            "will be SCE, AnnData, FlatFile, Seurat and HTAN. ")
        format <- c("SCE", "AnnData", "FlatFile", "HTAN", "Seurat")
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
        exportSCEtoFlatFile(inSCE, outputDir = fn, prefix=samplename)
    }

    if ("AnnData" %in% format) {
        ## Export to Python AnnData
        fp <- file.path(directory, samplename, "Python")
        dir.create(fp, showWarnings = TRUE, recursive = TRUE)
        fn <- file.path(fp, type)
        exportSCEtoAnnData(inSCE, outputDir=fn, compression='gzip', prefix=samplename)
    }

    if ("Seurat" %in% format) {
        ## Export to Seurat object
        fp <- file.path(directory, samplename, "Seurat")
        dir.create(fp, showWarnings = TRUE, recursive = TRUE)
        prefix <- paste0(samplename , paste0("_", type, ".rds"))
        exportSCEToSeurat(inSCE, prefix = prefix, outputDir = fp, overwrite = TRUE)
    }
}


#' Generate HTAN manifest file for droplet and cell count data
#' @param dropletSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' droplet count matrix data
#' @param cellSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' cell count matrix data
#' @param samplename The sample name of the \link[SingleCellExperiment]{SingleCellExperiment} objects
#' @param dir The output directory of the SCTK QC pipeline.
#' @param dataType Type of the input data. It can be one of "Droplet", "Cell" or "Both".
#' @param HTAN Whether generates manifest file including HTAN specific ID (HTAN Biospecimen ID, 
#' HTAN parent file ID and HTAN patient ID). Default is TRUE. 
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object which combines all
#' objects in sceList. The colData is merged.
#' @export

generateMeta <- function(dropletSCE = NULL,
                          cellSCE = NULL,
                          samplename,
                          dir,
                          HTAN=TRUE,
                          dataType = c("Droplet", "Cell", "Both")) {
  level3List <- list()
  level4List <- list()
  dataType = match.arg(dataType)

  directory <- file.path(basename(dir), samplename)
  filterDir <- file.path(directory, 'FlatFile', 'Cells')
  rawDir <- file.path(directory, 'FlatFile', 'Droplets')
  pkgVersion <- Biobase::package.version("singleCellTK")

  WorkFlowData = c(
    WorkFlow = 'Other',
    WorkFlowVer = paste('singleCellTK', pkgVersion, sep=':'),
    ParRaw = 'Ran perCellQC, EmptyDrops and barcodeRankDrops using singleCellTK',
    ParFiltered = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    PardecontX = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    PardecontX_bg = 'Ran perCellQC, doublet detection and decontX(with background matrix) using singleCellTK',
    ParSoupX = 'Ran perCellQC, doublet detection and SoupX using singleCellTK',
    ParSoupX_bg = 'Ran perCellQC, doublet detection and SoupX(with background matrix) using singleCellTK',
    ColData = 'Ran perCellQC, doublet detection and decontX using singleCellTK.',
    DecontXUMAP = 'UMAP dimension reduction generated by decontX.',
    DecontXUMAPBP = 'UMAP dimension reduction generated by decontX(with background matrix).',
    SoupXUMAP = 'UMAP dimension reduction generated by SoupX.',
    SoupXUMAPBP = 'UMAP dimension reduction generated by SoupX(with background matrix).',
    ScrubletTSNE = 'tSNE dimension reduction generated by Scrublet.',
    ScrubletUMAP = 'UMAP dimension reduction generated by Scrublet.'
  )
  


  ### calculate summary stats
  if (dataType == "Droplet" | dataType == "Both") {
    droplet_stat = c(CellNum = ncol(dropletSCE),
              MedianReads = stats::median(colData(dropletSCE)$sum),
              MedianGenes = stats::median(colData(dropletSCE)$detected),
              DataType = 'Droplet Matrix',
              FileName = file.path(rawDir, 'assays', paste0(samplename,'_counts.mtx.gz')))
  }

  if (dataType == "Cell" | dataType == "Both") {
    if ("decontXcounts" %in% SummarizedExperiment::assayNames(cellSCE)) {
      decontX_stat = c(CellNum = ncol(cellSCE),
                       MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts'))),
                       MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts') > 0)),
                       DataType = 'Decontaminated cell matrix return returned by runDecontX',
                       FileName = file.path(filterDir, 'assays', paste0(samplename,'_decontXcounts.mtx.gz')))      
    }

    if ("decontXcounts_bg" %in% SummarizedExperiment::assayNames(cellSCE)) {
      decontX_bg_stat = c(CellNum = ncol(cellSCE),
                MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts_bg'))),
                MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts_bg') > 0)),
                DataType = 'Decontaminated cell matrix returned by runDecontX , which is run with background count matrix',
                FileName = file.path(filterDir, 'assays', paste0(samplename,'_decontXcounts_bg.mtx.gz')))      
    }

    if ("SoupX" %in% SummarizedExperiment::assayNames(cellSCE)) {
      SoupX_stat = c(CellNum = ncol(cellSCE),
                       MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX'))),
                       MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX') > 0)),
                       DataType = 'Decontaminated cell matrix return returned by SoupX',
                       FileName = file.path(filterDir, 'assays', paste0(samplename,'_SoupX.mtx.gz')))      
    }

    if ("SoupX_bg" %in% SummarizedExperiment::assayNames(cellSCE)) {
      SoupX_bg_stat = c(CellNum = ncol(cellSCE),
                MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX_bg'))),
                MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX_bg') > 0)),
                DataType = 'Decontaminated cell matrix returned by SoupX , which is run with background count matrix',
                FileName = file.path(filterDir, 'assays', paste0(samplename,'_SoupX_bg.mtx.gz')))      
    }



    cell_stat = c(CellNum = ncol(cellSCE),
                   MedianReads = stats::median(colData(cellSCE)$sum),
                   MedianGenes = stats::median(colData(cellSCE)$detected),
                   DataType = 'Cell Matrix',
                   FileName = file.path(filterDir, 'assays', paste0(samplename,'_counts.mtx.gz')),
                   ColData = file.path(filterDir, paste0(samplename,'_cellData.txt.gz')),
                   DecontXUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP.txt.gz')),
                   DecontXUMAPBP = file.path(filterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP_bg.txt.gz')),
                   SoupXUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_SoupX_UMAP_all_cells.txt.gz')),
                   SoupXUMAPBP = file.path(filterDir, 'reducedDims', paste0(samplename,'_SoupX_bg_UMAP_all_cells.txt.gz')),
                   ScrubletTSNE = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_TSNE.txt.gz')),
                   ScrubletUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_UMAP.txt.gz')))
  }

  data <- list("Raw" = if (exists("droplet_stat")) {droplet_stat} else {NULL},
               "Filtered" = if (exists("cell_stat")) {cell_stat} else {NULL},
               "decontX" =  if (exists("decontX_stat")) {decontX_stat} else {NULL},
               "decontX_bg" =  if (exists("decontX_bg_stat")) {decontX_bg_stat} else {NULL},
               "SoupX" =  if (exists("SoupX_stat")) {SoupX_stat} else {NULL},
               "SoupX_bg" =  if (exists("SoupX_bg_stat")) {SoupX_bg_stat} else {NULL}
              )


  if (dataType %in% c("Cell", "Both")) {

    if (dataType == "Cell") {
      types <- c("Filtered")
    } else {
      types <- c("Filtered", "Raw")
    }
    
    l4Metrics <- c('ColData', 'ScrubletTSNE', 'ScrubletUMAP')
    if ('decontXcounts' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "decontX")
      l4Metrics <- c(l4Metrics, "DecontXUMAP")
    }

    if ('decontXcounts_bg' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "decontX_bg")
      l4Metrics <- c(l4Metrics, "DecontXUMAPBP")
    }

    if ('SoupX' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "SoupX")
      l4Metrics <- c(l4Metrics, "SoupXUMAP")
    }

    if ('SoupX_bg' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "SoupX_bg")
      l4Metrics <- c(l4Metrics, "SoupXUMAPBP")
    }    
  } else if (dataType == "Droplet") { 
    types <- c("Raw") 
  }


  for (type in types) {
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
      for (metric in l4Metrics) {
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

#' Generate HTAN manifest file for droplet and cell count data
#' @param dropletSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' droplet count matrix data
#' @param cellSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' cell count matrix data
#' @param samplename The sample name of the \link[SingleCellExperiment]{SingleCellExperiment} objects
#' @param htan_biospecimen_id The HTAN biospecimen id of the sample in \link[SingleCellExperiment]{SingleCellExperiment} object
#' @param dir The output directory of the SCTK QC pipeline.
#' @param dataType Type of the input data. It can be one of "Droplet", "Cell" or "Both".
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object which combines all
#' objects in sceList. The colData is merged.
#' @export
#' @importFrom rlang .data

generateHTANMeta <- function(dropletSCE = NULL,
                         cellSCE = NULL,
                         samplename,
                         htan_biospecimen_id,
                         dir,
                         dataType = c("Droplet", "Cell", "Both")) {
  level3List <- list()
  level4List <- list()
  dataType = match.arg(dataType)
  
  directory <- file.path(basename(dir), samplename)
  filterDir <- file.path(directory, 'FlatFile', 'Cells')
  rawDir <- file.path(directory, 'FlatFile', 'Droplets')
  
  absFilterDir <- file.path(dir, samplename, 'FlatFile', 'Cells')
  absRawDir <- file.path(dir, samplename, 'FlatFile', 'Droplets')
  
  pkgVersion <- Biobase::package.version("singleCellTK")
  htan_patient_id <- paste(stringr::str_split(htan_biospecimen_id, 
                                              '_', 
                                              simplify = TRUE)[seq(2)], 
                           collapse = '_')
  
  WorkFlowData = c(
    WorkFlow = 'Other',
    WorkFlowVer = paste('singleCellTK', pkgVersion, sep=':'),
    ParRaw = 'Ran perCellQC, EmptyDrops and barcodeRankDrops using singleCellTK',
    ParFiltered = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    PardecontX = 'Ran perCellQC, doublet detection and decontX using singleCellTK',
    PardecontX_bg = 'Ran perCellQC, doublet detection and decontX(with background matrix) using singleCellTK',
    ParSoupX = 'Ran perCellQC, doublet detection and SoupX using singleCellTK',
    ParSoupX_bg = 'Ran perCellQC, doublet detection and SoupX(with background matrix) using singleCellTK',
    ColData = 'Ran perCellQC, doublet detection and decontX using singleCellTK. Detailed paramters please see: ',
    DecontXUMAP = 'UMAP dimension reduction generated by decontX. Detailed paramters please see: ',
    DecontXUMAPBP = 'UMAP dimension reduction generated by decontX(with background matrix). Detailed paramters please see: ',
    SoupXUMAP = 'UMAP dimension reduction generated by SoupX. Detailed paramters please see: ',
    SoupXUMAPBP = 'UMAP dimension reduction generated by SoupX(with background matrix). Detailed paramters please see: ',
    ScrubletTSNE = 'tSNE dimension reduction generated by Scrublet. Detailed paramters please see: ',
    ScrubletUMAP = 'UMAP dimension reduction generated by Scrublet. Detailed paramters please see: '
  )
  
  ### calculate summary stats
  if (dataType == "Droplet" | dataType == "Both") {
    droplet_stat = c(CellNum = ncol(dropletSCE),
                     MedianReads = stats::median(colData(dropletSCE)$sum),
                     MedianGenes = stats::median(colData(dropletSCE)$detected),
                     DataType = 'Droplet Matrix',
                     FileName = file.path(rawDir, 'assays', paste0(samplename,'_counts.mtx.gz')),
                     AbsFileName = file.path(absRawDir, 'assays', paste0(samplename,'_counts.mtx.gz')))
  }
  
  if (dataType == "Cell" | dataType == "Both") {
    if ("decontXcounts" %in% SummarizedExperiment::assayNames(cellSCE)) {
      decontX_stat = c(CellNum = ncol(cellSCE),
                       MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts'))),
                       MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts') > 0)),
                       DataType = 'Decontaminated cell matrix return returned by runDecontX',
                       FileName = file.path(filterDir, 'assays', paste0(samplename,'_decontXcounts.mtx.gz')),
                       AbsFileName = file.path(absFilterDir, 'assays', paste0(samplename,'_decontXcounts.mtx.gz')))      
    }

    if ("decontXcounts_bg" %in% SummarizedExperiment::assayNames(cellSCE)) {
      decontX_bg_stat = c(CellNum = ncol(cellSCE),
                MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts_bg'))),
                MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'decontXcounts_bg') > 0)),
                DataType = 'Decontaminated cell matrix returned by runDecontX , which is run with background count matrix',
                FileName = file.path(filterDir, 'assays', paste0(samplename,'_decontXcounts_bg.mtx.gz')),
                AbsFileName = file.path(absFilterDir, 'assays', paste0(samplename,'_decontXcounts_bg.mtx.gz')))      
    }

    if ("SoupX" %in% SummarizedExperiment::assayNames(cellSCE)) {
      SoupX_stat = c(CellNum = ncol(cellSCE),
                       MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX'))),
                       MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX') > 0)),
                       DataType = 'Decontaminated cell matrix return returned by SoupX',
                       FileName = file.path(filterDir, 'assays', paste0(samplename,'_SoupX.mtx.gz')),
                       AbsFileName = file.path(absFilterDir, 'assays', paste0(samplename,'_SoupX.mtx.gz')))      
    }

    if ("SoupX_bg" %in% SummarizedExperiment::assayNames(cellSCE)) {
      SoupX_bg_stat = c(CellNum = ncol(cellSCE),
                MedianReads = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX_bg'))),
                MedianGenes = stats::median(Matrix::colSums(assay(cellSCE, 'SoupX_bg') > 0)),
                DataType = 'Decontaminated cell matrix returned by SoupX , which is run with background count matrix',
                FileName = file.path(filterDir, 'assays', paste0(samplename,'_SoupX_bg.mtx.gz')),
                AbsFileName = file.path(absFilterDir, 'assays', paste0(samplename,'_SoupX_bg.mtx.gz')))      
    }

    cell_stat = c(CellNum = ncol(cellSCE),
                  MedianReads = stats::median(colData(cellSCE)$sum),
                  MedianGenes = stats::median(colData(cellSCE)$detected),
                  DataType = 'Cell Matrix',
                  FileName = file.path(filterDir, 'assays', paste0(samplename,'_counts.mtx.gz')),
                  AbsFileName = file.path(absFilterDir, 'assays', paste0(samplename,'_counts.mtx.gz')),
                  ColData = file.path(filterDir, paste0(samplename,'_cellData.txt.gz')),
                  AbsColData = file.path(absFilterDir, paste0(samplename,'_cellData.txt.gz')),
                  DecontXUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP.txt.gz')),
                  AbsDecontXUMAP = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP.txt.gz')),
                  DecontXUMAPBP = file.path(filterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP_bg.txt.gz')),                   
                  AbsDecontXUMAPBP = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_decontX_UMAP_bg.txt.gz')),
                  SoupXUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_SoupX_UMAP_all_cells.txt.gz')),
                  AbsSoupXUMAP = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_SoupX_UMAP_all_cells.txt.gz')),
                  SoupXUMAPBP = file.path(filterDir, 'reducedDims', paste0(samplename,'_SoupX_bg_UMAP_all_cells.txt.gz')),                   
                  AbsSoupXUMAPBP = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_SoupX_bg_UMAP_all_cells.txt.gz')),
                  ScrubletTSNE = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_TSNE.txt.gz')),
                  AbsScrubletTSNE = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_scrublet_TSNE.txt.gz')),
                  ScrubletUMAP = file.path(filterDir, 'reducedDims', paste0(samplename,'_scrublet_UMAP.txt.gz')),
                  AbsScrubletUMAP = file.path(absFilterDir, 'reducedDims', paste0(samplename,'_scrublet_UMAP.txt.gz')))
  }
  
  data <- list("Raw" = if (exists("droplet_stat")) {droplet_stat} else {NULL},
               "Filtered" = if (exists("cell_stat")) {cell_stat} else {NULL},
               "decontX" =  if (exists("decontX_stat")) {decontX_stat} else {NULL},
               "decontX_bg" =  if (exists("decontX_bg_stat")) {decontX_bg_stat} else {NULL},
               "SoupX" =  if (exists("SoupX_stat")) {SoupX_stat} else {NULL},
               "SoupX_bg" =  if (exists("SoupX_bg_stat")) {SoupX_bg_stat} else {NULL}
              )


  if (dataType %in% c("Cell", "Both")) {

    if (dataType == "Cell") {
      types <- c("Filtered")
    } else {
      types <- c("Filtered", "Raw")
    }
    
    l4Metrics <- c('ColData', 'ScrubletTSNE', 'ScrubletUMAP')
    if ('decontXcounts' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "decontX")
      l4Metrics <- c(l4Metrics, "DecontXUMAP")
    }

    if ('decontXcounts_bg' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "decontX_bg")
      l4Metrics <- c(l4Metrics, "DecontXUMAPBP")
    }

    if ('SoupX' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "SoupX")
      l4Metrics <- c(l4Metrics, "SoupXUMAP")
    }

    if ('SoupX_bg' %in% SummarizedExperiment::assayNames(cellSCE)) {
      types <- c(types, "SoupX_bg")
      l4Metrics <- c(l4Metrics, "SoupXUMAPBP")
    }    
  } else if (dataType == "Droplet") { 
    types <- c("Raw") 
  }
  
  for (type in types) {
    level3List[[type]] <- data.frame(
      'Component' = 'ScRNA-seqLevel3',
      'Filename' = data[[type]]['FileName'],
      'FileLoc' = data[[type]]['AbsFileName'],
      'File Format' = 'mtx',
      'HTAN Biospecimen ID' = htan_biospecimen_id, #DATA_TYPE = data[[type]]['DataType'],
      'Data Category' = 'Gene Expression', 'Matrix Type' = 'Raw Counts',
      'Cell Median Number Reads' = data[[type]]['MedianReads'],
      'Cell Median Number Genes' = data[[type]]['MedianGenes'],
      'Cell Total' = data[[type]]['CellNum'],
      'scRNAseq Workflow Type' = WorkFlowData['WorkFlow'],
      'scRNAseq Workflow Parameters Description' = WorkFlowData[paste0('Par', type)],
      'Workflow Link' = paste0('http://sctk.camplab.net/v', pkgVersion, '/index.html'), 
      'Workflow Version' = WorkFlowData['WorkFlowVer'],
      'Workflow Start Datetime' = '',
      'Workflow End Datetime' = '',
      stringsAsFactors = FALSE, check.names=FALSE)
    
    
    level3List[[type]] <- cbind(level3List[[type]],
                                'HTAN Patient ID' = htan_patient_id, 'HTAN Parent Data File ID' = '',
                                'HTAN Data File ID' = '', 'Linked Matrices' = '')
  
    
    if (type == 'Filtered') {
      for (metric in l4Metrics) {
        level4List[[metric]] <- data.frame(
          'Component' = 'ScRNA-seqLevel4', "Filename" = data[[type]][metric],
          'FileLoc' = data[[type]][paste0('Abs', metric)],
          'File Format' = 'txt', 
          'HTAN Biospecimen ID' = htan_biospecimen_id, 
          'scRNAseq Workflow Type' = 'Other', #WorkFlowData[metric]
          'scRNAseq Workflow Parameters Description' = paste0(WorkFlowData[metric], file.path(directory, paste0(samplename, '_QCParameters.yaml'))),
          'Workflow Version' = WorkFlowData['WorkFlowVer'],
          'Workflow Link' = paste0('http://sctk.camplab.net/v', pkgVersion, '/index.html'),
          'Workflow Start Datetime' = '',
          'Workflow End Datetime' = '',
          stringsAsFactors = FALSE, check.names=FALSE)

        level4List[[metric]] <- cbind(level4List[[metric]], 'HTAN Patient ID' = htan_patient_id,
                                      'HTAN Data File ID' = '','HTAN Parent Data File ID' = '')
      
      }
    }
  }
  

  level3Meta <- do.call(base::rbind, level3List)
  level4Meta <- do.call(base::rbind, level4List)

  level3Meta$checkSumInt <- vapply(level3Meta$FileLoc, function(f) {
    hash <- tools::md5sum(f)
    .convert_hex_to_int(hash)
  }, character(1))
  
  level3Meta$`HTAN Data File ID` <- paste(level3Meta$`HTAN Patient ID`, level3Meta$checkSumInt, sep="_")
  linkMat <- dplyr::group_by(level3Meta, .data$`HTAN Biospecimen ID`) %>% dplyr::summarise(LinkedMat = paste(.data$`HTAN Data File ID`, collapse=","))
  level3Meta$`Linked Matrices` <- unlist(linkMat[match(level3Meta$`HTAN Biospecimen ID`, linkMat$`HTAN Biospecimen ID`), 'LinkedMat'])
  
  ParentFile <- level3Meta[grep("_counts.mtx.gz", unlist(level3Meta['FileLoc'])), c('HTAN Biospecimen ID', 'HTAN Data File ID')]
  level4Meta$`HTAN Parent Data File ID` <- unlist(ParentFile[match(level4Meta$`HTAN Biospecimen ID`, ParentFile$`HTAN Biospecimen ID`), 
                                                             'HTAN Data File ID'])
  
  level4Meta$checkSumInt <- vapply(level4Meta$FileLoc, function(f) {
    hash <- tools::md5sum(f)
    .convert_hex_to_int(hash)
  }, character(1))
  level4Meta$`HTAN Data File ID` <- paste(level4Meta$`HTAN Patient ID`, level4Meta$checkSumInt, sep="_")
  
  level3Meta <- level3Meta[,c("HTAN Biospecimen ID",
                              "Component",
                              "Filename",
                              "File Format",
                              "HTAN Parent Data File ID",
                              "HTAN Data File ID",
                              "Data Category",
                              "Matrix Type",
                              "Linked Matrices",
                              "Cell Median Number Reads",
                              "Cell Median Number Genes",
                              "Cell Total",
                              "scRNAseq Workflow Type",
                              "scRNAseq Workflow Parameters Description",
                              "Workflow Link",
                              "Workflow Version",
                              "Workflow Start Datetime",
                              "Workflow End Datetime")]
  
  level4Meta <- level4Meta[,c("HTAN Biospecimen ID",
                              "Component",
                              "Filename",
                              "File Format",
                              "HTAN Parent Data File ID",
                              "HTAN Data File ID",
                              "scRNAseq Workflow Type",
                              "scRNAseq Workflow Parameters Description",
                              "Workflow Version",
                              "Workflow Link",
                              "Workflow Start Datetime",
                              "Workflow End Datetime")]
  

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
#' @return If \code{writeYAML} TRUE, a yaml object will be generated. If FALSE, character object.
#' @export
getSceParams <- function(inSCE,
                         skip = c("scrublet", "runDecontX","runBarcodeRanksMetaOutput"),
                         ignore = c("algorithms", "estimates","contamination",
                                    "z","sample","rank","BPPARAM","batch","geneSetCollection",
                                    "barcodeArgs"),
                         directory = './',
                         samplename = '',
                         writeYAML = TRUE) {

  meta <- S4Vectors::metadata(inSCE)
  algos <- names(meta)[!names(meta) %in% skip]
  outputs <- '---'
  parList <- list()
  dir <- file.path(directory, samplename)

  for (algo in algos) {
    params <- meta[[algo]][[1]]
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

    if (preproc == "Alevin") {
        cellSCE <- importAlevin(alevinDir = path, sampleName = samplename, class = "Matrix", delayedArray=FALSE)
        return(list(dropletSCE, cellSCE))
    }
    ## preproc is not one of the method above. Stop the pipeline.
    stop(paste0("'", preproc, "' not supported."))
}

