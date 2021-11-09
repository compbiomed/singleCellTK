
#' Imports samples from different sources and compiles them into a list of SCE objects
#' @param allImportEntries object containing the sources and parameters of all the samples being imported (from the UI)
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link{DelayedArray} object or not. Default \code{FALSE}.
#' @return A list of \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the droplet or cell data or both,depending on the dataType that users provided.
#' @export
importMultipleSources <- function(allImportEntries, delayedArray = FALSE) {
  sceObjs <- list()
  for (entry in allImportEntries$samples) {
    if (entry$type == "cellRanger2") {
      if (is.null(entry$params$cellRangerDirs)) {
        newSce <- importCellRangerV2Sample(
          dataDir = entry$params$dataDir,
          sampleName = entry$params$sampleName,
          delayedArray = delayedArray
        )
      } else {
        newSce <- importCellRangerV2(
          cellRangerDirs = entry$params$cellRangerDirs,
          sampleDirs = entry$params$sampleDirs,
          sampleNames = entry$params$sampleNames,
          delayedArray = delayedArray,
          reference = entry$params$reference
        )
      }

    } else if (entry$type == "cellRanger3") {
      if (is.null(entry$params$cellRangerDirs)) {
        newSce <- importCellRangerV3Sample(
          dataDir = entry$params$dataDir,
          sampleName = entry$params$sampleName,
          delayedArray = delayedArray
        )
      } else {
        newSce <- importCellRangerV3(
          cellRangerDirs = entry$params$cellRangerDirs,
          sampleDirs = entry$params$sampleDirs,
          sampleNames = entry$params$sampleNames,
          delayedArray = delayedArray
        )
      }
    } else if (entry$type == "starSolo") {
      newSce <- importSTARsolo(
        STARsoloDirs = entry$params$STARsoloDirs,
        samples = entry$params$amples,
        delayedArray = delayedArray
      )
    } else if (entry$type == "busTools") {
      newSce <- importBUStools(
        BUStoolsDirs = entry$params$BUStoolsDirs,
        samples = entry$params$samples,
        delayedArray = delayedArray
      )
    } else if (entry$type == "seqc") {
      newSce <- importSEQC(
        seqcDirs = entry$params$seqcDirs,
        samples = entry$params$samples,
        prefix = entry$params$prefix,
        delayedArray = delayedArray
      )
    } else if (entry$type == "optimus") {
      newSce <- importOptimus(
        OptimusDirs = entry$params$OptimusDirs,
        samples = entry$params$samples,
        delayedArray = delayedArray
      )
    } else if (entry$type == "files") {
      newSce <- importFromFiles(assayFile = entry$params$assayFile,
                                annotFile = entry$params$annotFile,
                                featureFile = entry$params$featureFile,
                                assayName = entry$params$assayName,
                                delayedArray = delayedArray)
    } else if (entry$type == "example") {
      newSce <- importExampleData(dataset = entry$params$dataset,
                                  delayedArray = delayedArray)
    } else if (entry$type == "rds") {
      importedrds <- readRDS(entry$params$rdsFile)
      if (base::inherits(importedrds, "SummarizedExperiment")) {
        newSce <- importedrds
      } else if (base::inherits(importedrds, "Seurat")) {
        newSce <- convertSeuratToSCE(importedrds)
      } else {
        stop("The '.rds' file should contain a 'SingleCellExperiment' or 'Seurat' object.")
      }

      for(assay in SummarizedExperiment::assayNames(newSce)) {
        if(!base::inherits(SummarizedExperiment::assay(newSce, assay), "dgCMatrix") && !isTRUE(delayedArray)) {
          SummarizedExperiment::assay(newSce, assay) <- .convertToMatrix(SummarizedExperiment::assay(newSce, assay))
        }

        if(!base::inherits(SummarizedExperiment::assay(newSce, assay), "DelayedArray") && isTRUE(delayedArray)) {
          SummarizedExperiment::assay(newSce, assay) <- DelayedArray::DelayedArray(SummarizedExperiment::assay(newSce, assay))
        }
      }

    }
    
    # Begin Set Tags
    if(entry$type %in% c("cellRanger2", "cellRanger3", "starSolo", "busTools", "seqc", "optimus", "example")){
      newSce <- expSetDataTag(
        inSCE = newSce,
        assayType = "raw",
        assays = SummarizedExperiment::assayNames(newSce))
    }
    else if(entry$type %in% c("rds", "files")){
      # Check if tags already stored in uploaded rds/files
      if(is.null(S4Vectors::metadata(newSce)$assayType)){
        try({
          counts(newSce)
          newSce <- expSetDataTag(
            inSCE = newSce,
            assayType = "raw",
            assays = "counts")
        }, silent = TRUE)
        
        try({
          logcounts(newSce)
          newSce <- expSetDataTag(
            inSCE = newSce,
            assayType = "transformed",
            assays = "logcounts")
        }, silent = TRUE)
        
        try({
          normcounts(newSce)
          newSce <- expSetDataTag(
            inSCE = newSce,
            assayType = "transformed",
            assays = "normcounts")
        }, silent = TRUE)
        
        try({
          celda::decontXcounts(newSce)
          newSce <- expSetDataTag(
            inSCE = newSce,
            assayType = "raw",
            assays = "decontXcounts")
        }, silent = TRUE)
        
        untaggedAssays <- SummarizedExperiment::assayNames(newSce)
        untaggedAssays <- untaggedAssays[! untaggedAssays %in% c('counts', 'logcounts', 'normcounts', 'decontX')]
        
        newSce <- expSetDataTag(
          inSCE = newSce,
          assayType = "uncategorized",
          assays = untaggedAssays) 
      }
      # End Set Tags
    }
    
    sceObjs = c(sceObjs, list(newSce))
  }
  
  return(combineSCE(sceList = sceObjs,
                    by.r = Reduce(base::intersect, lapply(sceObjs, function(x) { colnames(rowData(x))})),
                    by.c = Reduce(base::intersect, lapply(sceObjs, function(x) { colnames(colData(x))})),
                    combined = TRUE)
  )
}

