
#' Imports samples from different sources and compiles them into a list of SCE objects
#' @param allImportEntries object containing the sources and parameters of all the samples being imported (from the UI)
#' @return A list of \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the droplet or cell data or both,depending on the dataType that users provided.
#' @export
importMultipleSources <- function(allImportEntries) {
  sceObjs <- list()
  for (entry in allImportEntries$samples) {
    if (entry$type == "cellRanger2") {
      if (is.null(entry$params$cellRangerDirs)) {
        newSce <- importCellRangerV2Sample(
          dataDir = entry$params$dataDir,
          sampleName = entry$params$sampleName,
        )
      } else {
        newSce <- importCellRangerV2(
          cellRangerDirs = entry$params$cellRangerDirs,
          sampleDirs = entry$params$sampleDirs,
          sampleNames = entry$params$sampleNames,
        )
      }
      
    } else if (entry$type == "cellRanger3") {
      if (is.null(entry$params$cellRangerDirs)) {
        newSce <- importCellRangerV3Sample(
          dataDir = entry$params$dataDir,
          sampleName = entry$params$sampleName,
        )
      } else {
        newSce <- importCellRangerV3(
          cellRangerDirs = entry$params$cellRangerDirs,
          sampleDirs = entry$params$sampleDirs,
          sampleNames = entry$params$sampleNames,
        )
      }
    } else if (entry$type == "starSolo") {
      newSce <- importSTARsolo(
        STARsoloDirs = entry$params$STARsoloDirs,
        samples = entry$params$amples
      )
    } else if (entry$type == "busTools") {
      newSce <- importBUStools(
        BUStoolsDirs = entry$params$BUStoolsDirs,
        samples = entry$params$samples,
      )
    } else if (entry$type == "seqc") {
      newSce <- importSEQC(
        seqcDirs = entry$params$seqcDirs,
        samples = entry$params$samples,
        prefix = entry$params$prefix,
      )
    } else if (entry$type == "optimus") {
      newSce <- importOptimus(
        OptimusDirs = entry$params$OptimusDirs,
        samples = entry$params$samples
      )
    } else if (entry$type == "files") {
      newSce <- importFromFiles(assayFile = entry$params$assayFile,
                                annotFile = entry$params$annotFile,
                                featureFile = entry$params$featureFile,
                                assayName = entry$params$assayName)
    } else if (entry$type == "example") {
      newSce <- importExampleData(dataset = entry$params$dataset)
    } else if (entry$type == "rds") {
      importedrds <- readRDS(entry$params$rdsFile)
      if (base::inherits(importedrds, "SummarizedExperiment")) {
        newSce <- importedrds
      } else if (base::inherits(importedrds, "Seurat")) {
        newSce <- convertSeuratToSCE(importedrds)
      } else {
        stop("The '.rds' file should contain a 'SingleCellExperiment' or 'Seurat' object.")
      }
    }
    sceObjs = c(sceObjs, list(newSce))
  }
  
  return(combineSCE(sceList = sceObjs,
                    by.r = NULL,
                    by.c = Reduce(intersect, lapply(sceObjs, function(x) { colnames(colData(x))})),
                    combined = TRUE)
  )
}

