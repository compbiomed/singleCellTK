#' @title sceTxtExport
#' @description exports SingleCellExperiment object to txt.gz file(s).
#' @param sce SingleCellExperiment object to be exported.
#' @param outputDir Name of the directory to store the exported txt.gz file(s).
#' @examples
#' sceTxtExport(sce, outputDir)
#' @export txt.gz file(s)
#'
library(SingleCellExperiment)
library(SummarizedExperiment)

sceTxtExport <- function(sce, outputDir) {
  # creating folders
  assays_folder <- paste0(outputDir, "/assays")
  
  if (dir.exists(outputDir)) {
    print("Warning: Output folder already exists. Overwriting.")
  } else {
    dir.create(outputDir)
  }
  
  if (!dir.exists(assays_folder)) {
    dir.create(assays_folder)
  }
  
  if (length(reducedDimNames(sce)) > 0) {
    reduce_dimensions_folder <- paste0(outputDir, "/reducedDimNames")
    
    if (!dir.exists(reduce_dimensions_folder)) {
      dir.create(reduce_dimensions_folder)
    }
  }
  
  if (length(altExpNames(sce)) > 0) {
    altExp_folder <- paste0(outputDir, "/altExp")
    if (!dir.exists(altExp_folder)) {
      dir.create(altExp_folder)
    }
  }
  
  
  
  if (length(names(metadata(sce))) > 0) {
    metadata_folder <- paste0(outputDir, "/metadata")
    if (!dir.exists(metadata_folder)) {
      dir.create(metadata_folder)
    }
  }
  
  # function to write txt gz files
  writeGzFile <- function(data, name) {
    data <- as.data.frame(data)
    gz <- gzfile(paste0(name, ".txt.gz"), open = "w")
    write.table(as.data.frame(data), gz, row.names = T)
    close(con = gz, type = "w")
  }
  
  # write assays data
  assay_names <- names(assays(sce))
  for (i in assay_names) {
    data <- assays(sce)[[i]]
    filename <- paste0(assays_folder, "/", i)
    print(paste("Writing", filename, "..."))
    writeGzFile(data, filename)
  }
  
  # write altExpNames
  if (length(altExpNames(sce)) > 0) {
    altExp_names <- altExpNames(sce)
    for (i in altExp_names) {
      data <- altExp(sce)[[i]]
      filename <- paste0(altExp_folder, "/", i)
      print(paste("Writing", filename, "..."))
      writeGzFile(data, filename)
    }
  }
  
  # write colData and rowData
  
  # write colData
  data <- colData(sce)
  if (length(data) > 0) {
    filename <- paste0(outputDir, "/", "colData")
    print(paste("Writing", filename, "..."))
    writeGzFile(data, name = filename)
  }
  
  # write rowData
  
  data <- rowData(sce)
  if (length(colData(sce)) > 0) {
    filename <- paste0(outputDir, "/", "rowData")
    print(paste("Writing", filename, "..."))
    writeGzFile(data, filename)
  }
  
  
  
  # write reduced dimensions
  if (length(reducedDimNames(sce)) != 0) {
    redDim_names <- reducedDimNames(sce)
    for (i in redDim_names) {
      data <- reducedDim(sce, i, withDimnames = TRUE)
      filename <- paste0(reduce_dimensions_folder, "/", i)
      print(paste("Writing", filename, "..."))
      writeGzFile(data, filename)
    }
  }
  
  # write metaData
  
  if (length(names(metadata(sce))) > 0) {
    metadata_names <- names(metadata(sce))
    for (i in metadata_names) {
      data <- metadata(sce)[[i]]
      filename <- paste0(metadata_folder, "/", i)
      print(paste("Writing", filename, "..."))
      writeGzFile(data, filename)
    }
  }
  
  print("Done")
  
}

# testing

testFUN() {
  library(scRNAseq)
  
  data <- readRDS("data/pbmc3k_unfiltered.rds")
  sce <-
    SingleCellExperiment(assays = list(counts = data, logcounts = data))
  
  sce1 <- ReprocessedAllenData("tophat_counts")
  
  sceTxtExport(sce)
  sceTxtExport(sce1)
  
}