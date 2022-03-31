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
#' @param prefix Prefix of file names.
#' @return Generates text files containing data from \code{inSCE}.
#' @examples
#' data(sce_chcl, package = "scds")
#' \dontrun{
#' exportSCEtoFlatFile(sce_chcl, "sce_chcl")
#' }
#' @export
#' @importFrom SummarizedExperiment colData rowData
exportSCEtoFlatFile <- function(sce,
                                outputDir = "./",
                                overwrite = TRUE,
                                gzipped = TRUE,
                                prefix = 'SCE') {
  #path <- file.path(outputDir, prefix)
  path <- outputDir
  if (!file.exists(path)){
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }
  .writeAssays(sce, path, overwrite, gzipped, prefix)
  .writeColData(sce, path, overwrite, gzipped, prefix)
  .writeRowData(sce, path, overwrite, gzipped, prefix)
  .writeMetaData(sce, path, overwrite, prefix)
  .writeReducedDims(sce, path, overwrite, gzipped, prefix)
  .writeAltExps(sce, path, overwrite, gzipped, prefix)

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
  message(date(), " .. Writing '", filename)
  .checkOverwrite(filename, overwrite)
  data.table::fwrite(x = data, file = filename, nThread = 1, row.names = FALSE)
}


# write assays data
.writeAssays <- function(sce, path, overwrite, gzipped, sample) {
  if (length(SummarizedExperiment::assays(sce)) > 0) {
    assaysFolder <- file.path(path, "/assays")
    dir.create(assaysFolder, showWarnings = FALSE, recursive = TRUE)

    assayNames <- names(SummarizedExperiment::assays(sce))
    if(is.null(assayNames)) {
      assayNames <- paste0("assay", seq(SummarizedExperiment::assays(sce)))
    }
    for (i in seq_along(SummarizedExperiment::assays(sce))) {

      filename <- paste(sample, paste0(assayNames[i], ".mtx"), sep="_")
      assaypath <- file.path(assaysFolder, filename)

      .checkOverwrite(assaypath, overwrite)
      mat <- .convertToMatrix(SummarizedExperiment::assays(sce)[[i]])
      message(date(), " .. Writing assay '", assayNames[i], "' to ", assaypath)
      out <- Matrix::writeMM(mat, assaypath)

      if(isTRUE(gzipped)) {
        .checkOverwrite(paste0(assaypath, ".gz"), overwrite)
        message(date(), " .. Compressing into ", paste0(assaypath, ".gz"))
        R.utils::gzip(filename = assaypath, overwrite = overwrite)
      }
    }
  }
}



# write altExpNames
.writeAltExps <- function(sce, path, overwrite, gzipped, sample) {
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
      .writeAssays(sceAltExp, path = assaysFolder, overwrite = overwrite, gzipped = gzipped, sample)
      .writeColData(sceAltExp, altExpPath, overwrite, gzipped, sample)
      .writeRowData(sceAltExp, altExpPath, overwrite, gzipped, sample)
      .writeMetaData(sceAltExp, altExpPath, overwrite, sample)
    }
  }
}

# write colData
.writeColData <- function(sce, path, overwrite, gzipped, sample) {
  if(ncol(colData(sce)) > 0) {
    data <- SummarizedExperiment::colData(sce)
    colDataPath <- file.path(path, paste(sample, "cellData", sep="_"))
    .writeSCEFile(data, colDataPath, overwrite, gzipped)
  }
}


# write rowData
.writeRowData <- function(sce, path, overwrite, gzipped, sample) {
  if(ncol(rowData(sce)) > 0) {
    data <- SummarizedExperiment::rowData(sce)
    rowDataPath <-  file.path(path, paste(sample, "featureData", sep="_"))
    .writeSCEFile(data, rowDataPath, overwrite, gzipped)
  }
}


# write reduced dimensions
.writeReducedDims <- function(sce, path, overwrite, gzipped, sample) {
  if (length(SingleCellExperiment::reducedDimNames(sce)) > 0) {
    reducedDimsFolder <- file.path(path, "reducedDims")
    dir.create(reducedDimsFolder, showWarnings = FALSE, recursive = TRUE)

    if (length(reducedDimNames(sce)) > 0) {
      reducedDimNames <- SingleCellExperiment::reducedDimNames(sce)
      for (i in reducedDimNames) {
        data <- SingleCellExperiment::reducedDim(sce, i, withDimnames = TRUE)
        reducedDimNamePath <- file.path(reducedDimsFolder, paste(sample, i, sep="_"))
        message(date(), " .. Writing reducedDim '", i, "' to", reducedDimNamePath)
        .writeSCEFile(data, reducedDimNamePath, overwrite, gzipped)
      }
    }
  }
}

# write metaData
.writeMetaData <- function(sce, path, overwrite, sample) {
  if (length(S4Vectors::metadata(sce)) > 0) {
    metadataFolder <- file.path(path, "/metadata")
    dir.create(metadataFolder, showWarnings = FALSE, recursive = TRUE)

    filename <- file.path(metadataFolder, paste(sample, "metadata.rds", sep="_"))
    .checkOverwrite(filename, overwrite)
    message(date(), " .. Writing metadata to ", filename)
    saveRDS(object = S4Vectors::metadata(sce), file = filename)
  }
}

