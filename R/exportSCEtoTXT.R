#' @title Export a \link[SingleCellExperiment]{SingleCellExperiment} object to flat text files
#' @description Writes all assays, colData, rowData, reducedDims, and altExps objects in a
#' \link[SingleCellExperiment]{SingleCellExperiment} to text files.
#' The items in the 'metadata' slot remain stored in list and are saved in an RDS file.
#' @param sce \link[SingleCellExperiment]{SingleCellExperiment} object to be
#'  exported.
#' @param outputDir Name of the directory to store the exported file(s).
#' @param overwrite Boolean. Whether to overwrite the output files. Default
#'  \code{FALSE}.
#' @param outputType The desired type of the output files. Default \code{"txt"}
#'  which
#'  writes the \link[SingleCellExperiment]{SingleCellExperiment} object as
#'  tab delimited text files.
#' @param gzipped Boolean. \code{TRUE} if the output files are to be
#'  gzip compressed. \code{FALSE} otherwise. Default
#'  \code{TRUE} to save disk space.
#' @examples
#' data(sce_chcl, package = "scds")
#' exportSCEtoTXT(sce_chcl, "sce_chcl")
#'
#' @export
exportSCEtoTXT <- function(sce,
                     outputDir = "./",
                     outputType = "txt",
                     overwrite = FALSE,
                     gzipped = TRUE) {
  
  if (length(SummarizedExperiment::assays(sce)) > 0) {
    assaysFolder <- file.path(outputDir, "/assays")
    dir.create(assaysFolder, showWarnings = FALSE, recursive = TRUE)
    
    .writeAssays(sce, assaysFolder, overwrite, gzipped)
  }
  
  .writeColData(sce, outputDir, overwrite, gzipped)
  .writeRowData(sce, outputDir, overwrite, gzipped)
  
  if (length(SingleCellExperiment::reducedDimNames(sce)) > 0) {
    reducedDimsFolder <- file.path(outputDir, "/reducedDims")
    dir.create(reducedDimsFolder, showWarnings = FALSE, recursive = TRUE)
    .writeReducedDims(sce, reducedDimsFolder, overwrite, gzipped)
  }
  
  if (length(SingleCellExperiment::altExpNames(sce)) > 0) {
    altExpsFolder <- file.path(outputDir, "/altExps")
    dir.create(altExpsFolder, showWarnings = FALSE, recursive = TRUE)
    .writeAltExps(sce, altExpsFolder, overwrite, gzipped)
  }
  
}

.checkOverwrite <- function(path, overwrite) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(paste0(path, " already exists. Change 'outputDir' or set 'overwrite' to TRUE."))
  }
}

# function to write txt gz files
.writeSCEFile <- function(data, path, overwrite, gzipped) {
  data <- data.table::as.data.table(data, keep.rownames = TRUE, key=NULL, sorted=TRUE,
                                    value.name="value", na.rm=TRUE)
  if (isTRUE(gzipped)) {
    filename <- paste0(path, ".txt.gz")
  } else {
    filename <- paste0(path, ".txt")
  }
  
  .checkOverwrite(filename, overwrite)
  data.table::fwrite(data, file = filename)

}


# write assays data
.writeAssays <- function(sce, path, overwrite, gzipped) {
   assayNames <- names(SummarizedExperiment::assays(sce))
    for (i in assayNames) {
      data <- SummarizedExperiment::assays(sce)[[i]]
      assaypath <- file.path(path, i)
      message(date(), " Writing assay '", i, "'")
      
      match <- intersect(c("dgTMatrix", "dgCMatrix", "dgRMatrix", "matrix"), class(data))
      if (length(match) > 0) {
        filename <- paste0(assaypath, ".mtx")
        Matrix::writeMM(Matrix::Matrix(data, sparse = TRUE), filename)

    if (!("data.frame" %in% class(data))) {
      data <- data.table::as.data.table(Matrix::as.matrix(data), keep.rownames = TRUE)
      
    }

    .writeSCEFile(data, assaypath, overwrite, gzipped)
  }
    }
}


# write altExpNames
.writeAltExps <- function(sce, path, overwrite, gzipped) {
  altExpNames <- SingleCellExperiment::altExpNames(sce)
  for (i in altExpNames) {
    sceAltExp <- SingleCellExperiment::altExp(sce, i, withColData = FALSE)
    altExpPath <- file.path(path, i)
    message(date(), " Writing altExp '", i, "'")

    if (length(SummarizedExperiment::assays(sceAltExp)) > 0) {
      assaysFolder <- file.path(altExpPath, "/assays")
      dir.create(assaysFolder, showWarnings = FALSE, recursive = TRUE)
      .writeAssays(sceAltExp, assaysFolder, overwrite, gzipped)
    }

    .writeColData(sceAltExp, altExpPath, overwrite, gzipped)
    .writeRowData(sceAltExp, altExpPath, overwrite, gzipped)

    if (length(SingleCellExperiment::reducedDimNames(sceAltExp)) > 0) {
      reducedDimsFolder <- file.path(altExpPath, "/reducedDims")
      dir.create(reducedDimsFolder, showWarnings = FALSE, recursive = TRUE)
      .writeReducedDims(sceAltExp, reducedDimsFolder, overwrite, gzipped)
    }

    if (length(SingleCellExperiment::altExpNames(sceAltExp)) > 0) {
      altExpsFolder <- file.path(altExpPath, "/altExps")
      dir.create(altExpsFolder, showWarnings = FALSE, recursive = TRUE)
      .writeAltExps(sceAltExp, altExpsFolder, overwrite, gzipped)
    }

    if (length(names(S4Vectors::metadata(sceAltExp))) > 0) {
      metadataFolder <- file.path(altExpPath, "/metadata")
      dir.create(metadataFolder, showWarnings = FALSE, recursive = TRUE)
      .writeMetaData(sce, metadataFolder, overwrite)
    }
  }
}


# write colData
.writeColData <- function(sce, path, overwrite, gzipped) {
  data <- SummarizedExperiment::colData(sce)
  colDataPath <- file.path(path, "colData")
  message(date(), " Writing colData")
  .writeSCEFile(data, colDataPath, overwrite, gzipped)
}


# write rowData
.writeRowData <- function(sce, path, overwrite, gzipped) {
  data <- SummarizedExperiment::rowData(sce)
  rowDataPath <- file.path(path, "rowData")
  message(date(), " Writing rowData")
  .writeSCEFile(data, rowDataPath, overwrite, gzipped)
}


# write reduced dimensions
.writeReducedDims <- function(sce, path, overwrite, gzipped) {
  if (length(reducedDimNames(sce)) > 0) {
    reducedDimNames <- SingleCellExperiment::reducedDimNames(sce)
    for (i in reducedDimNames) {
      data <- SingleCellExperiment::reducedDim(sce, i, withDimnames = TRUE)
      reducedDimNamePath <- file.path(path, i)
      message(date(), " Writing reducedDim '", i, "'")
      .writeSCEFile(data, reducedDimNamePath, overwrite, gzipped)
    }
  }
}


# write metaData
.writeMetaData <- function(sce, outputDir, overwrite) {
  if (length(S4Vectors::metadata(sce)) > 0) {
    metadataFolder <- file.path(outputDir, "/metadata")
    dir.create(metadataFolder, showWarnings = FALSE, recursive = TRUE)
    saveRDS(object = S4Vectors::metadata(sce), file = file.path(metadataFolder, "metadata.rds"))
  }
}


