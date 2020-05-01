#' @title Export a \link[SingleCellExperiment]{SingleCellExperiment} object to flat text files
#' @description Writes all assays, colData, rowData, reducedDims, and altExps objects in a
#' \link[SingleCellExperiment]{SingleCellExperiment} to text files.
#' The items in the 'metadata' slot remain stored in list and are saved in an RDS file.
#' @param sce \link[SingleCellExperiment]{SingleCellExperiment} object to be
#'  exported.
#' @param outputDir Name of the directory to store the exported file(s).
#' @param overwrite Boolean. Whether to overwrite the output files. Default
#'  \code{TRUE}.
#' @param gzipped Boolean. \code{TRUE} if the output files are to be
#'  gzip compressed. \code{FALSE} otherwise. Default
#'  \code{TRUE}.
#' @examples
#' data(sce_chcl, package = "scds")
#' exportSCEtoFlatFile(sce_chcl, "sce_chcl")
#'
#' @export
exportSCEtoFlatFile <- function(sce,
                     outputDir = "./",
                     overwrite = TRUE,
                     gzipped = TRUE) {
  
  .writeAssays(sce, outputDir, overwrite, gzipped)
  .writeColData(sce, outputDir, overwrite, gzipped)  
  .writeRowData(sce, outputDir, overwrite, gzipped)  
  .writeMetaData(sce, outputDir, overwrite)
  .writeReducedDims(sce, outputDir, overwrite, gzipped)
  .writeAltExps(sce, outputDir, overwrite, gzipped)

}

.checkOverwrite <- function(path, overwrite) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(paste0(path, " already exists. Change 'outputDir' or set 'overwrite' to TRUE."))
  }
}

# function to write txt gz files
.writeSCEFile <- function(data, path, overwrite, gzipped) {

  if(is.null(rownames(data))) {
    data <- data.frame(ID = seq(nrow(data)), data)
  } else {
    data <- data.frame(ID = rownames(data), data)
  }
  
  if (isTRUE(gzipped)) {
    filename <- paste0(path, ".txt.gz")
  } else {
    filename <- paste0(path, ".txt")
  }
  
  .checkOverwrite(filename, overwrite)
  data.table::fwrite(x = data, file = filename, nThread = 1, row.names = FALSE)
}


# write assays data
.writeAssays <- function(sce, path, overwrite, gzipped) {
  if (length(SummarizedExperiment::assays(sce)) > 0) {
    assaysFolder <- file.path(path, "/assays")
    dir.create(assaysFolder, showWarnings = FALSE, recursive = TRUE)
    
    assayNames <- names(SummarizedExperiment::assays(sce))
    if(is.null(assayNames)) {
      assayNames <- paste0("assay", seq(SummarizedExperiment::assays(sce)))
    }
    for (i in seq_along(SummarizedExperiment::assays(sce))) {
      message(date(), " .. Writing assay '", assayNames[i], "'")
      assaypath <- file.path(assaysFolder, paste0(assayNames[i], ".mtx"))
      
      .checkOverwrite(assaypath, overwrite)
      mat <- .convertToMatrix(SummarizedExperiment::assays(sce)[[i]])
      out <- Matrix::writeMM(mat, assaypath)
    
      if(isTRUE(gzipped)) {
        .checkOverwrite(paste0(assaypath, ".gz"), overwrite)
        R.utils::gzip(filename = assaypath, overwrite = overwrite)
      }
    }  
  }
}



# write altExpNames
.writeAltExps <- function(sce, path, overwrite, gzipped) {
  if (length(SingleCellExperiment::altExpNames(sce)) > 0) {
    altExpsFolder <- file.path(path, "/altExps")
    dir.create(altExpsFolder, showWarnings = FALSE, recursive = TRUE)
    
    altExpNames <- SingleCellExperiment::altExpNames(sce)
    for (i in altExpNames) {
      sceAltExp <- SingleCellExperiment::altExp(sce, i, withColData = FALSE)
      altExpPath <- file.path(path, i)
      message(date(), " .. Writing altExp '", i, "'")

      assaysFolder <- file.path(altExpPath, "/assays")
      dir.create(assaysFolder, showWarnings = FALSE, recursive = TRUE)
      .writeAssays(sceAltExp, path = assaysFolder, overwrite = overwrite, gzipped = gzipped)
      .writeColData(sceAltExp, altExpPath, overwrite, gzipped)
      .writeRowData(sceAltExp, altExpPath, overwrite, gzipped)
      .writeMetaData(sceAltExp, altExpPath, overwrite)
    }
  }
}

# write colData
.writeColData <- function(sce, path, overwrite, gzipped) {
  
  if(ncol(colData(sce)) > 0) {
    data <- SummarizedExperiment::colData(sce)
    colDataPath <- file.path(path, "colData")
    .writeSCEFile(data, colDataPath, overwrite, gzipped)
  } 
}


# write rowData
.writeRowData <- function(sce, path, overwrite, gzipped) {
  if(ncol(rowData(sce)) > 0) {
    data <- SummarizedExperiment::rowData(sce)
    rowDataPath <- file.path(path, "rowData")
    .writeSCEFile(data, rowDataPath, overwrite, gzipped)
  }
}


# write reduced dimensions
.writeReducedDims <- function(sce, path, overwrite, gzipped) {
  if (length(SingleCellExperiment::reducedDimNames(sce)) > 0) {
    reducedDimsFolder <- file.path(path, "/reducedDims")
    dir.create(reducedDimsFolder, showWarnings = FALSE, recursive = TRUE)
    
    if (length(reducedDimNames(sce)) > 0) {
      reducedDimNames <- SingleCellExperiment::reducedDimNames(sce)
      for (i in reducedDimNames) {
        message(date(), " .. Writing reducedDim '", i, "'")
        data <- SingleCellExperiment::reducedDim(sce, i, withDimnames = TRUE)
        reducedDimNamePath <- file.path(reducedDimsFolder, i)
        .writeSCEFile(data, reducedDimNamePath, overwrite, gzipped)
      }
    }
  }
}

# write metaData
.writeMetaData <- function(sce, path, overwrite) {
  if (length(S4Vectors::metadata(sce)) > 0) {
    metadataFolder <- file.path(path, "/metadata")
    dir.create(metadataFolder, showWarnings = FALSE, recursive = TRUE)
    
    filename <- file.path(metadataFolder, "metadata.rds")
    .checkOverwrite(filename, overwrite)
    saveRDS(object = S4Vectors::metadata(sce), file = filename)
  }
}


