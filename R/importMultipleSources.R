
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
    } else if (entry$type == "cellRanger3_files") {
      
      # get current tempDir by shiny
      mytempdir <- tempdir(check=TRUE)
      
      if(dir.exists(mytempdir)){
        # get uploaded filepaths
        matrixFilePath <- entry$params$assayFile
        barcodesFilePath <- entry$params$annotFile
        featuresFilePath <- entry$params$featureFile
        metricsFilePath <- entry$params$summaryFile
        
        # rename to original names
        file.rename(matrixFilePath, paste0(dirname(matrixFilePath), "/matrix.mtx.gz"))
        file.rename(barcodesFilePath, paste0(dirname(barcodesFilePath), "/barcodes.tsv.gz"))
        file.rename(featuresFilePath, paste0(dirname(featuresFilePath), "/features.tsv.gz"))
        if(!is.null(metricsFilePath)){
          file.rename(metricsFilePath, paste0(dirname(metricsFilePath), "/metrics_summary.csv"))
        }
        
        # create a sample folder
        dir.create(paste0(mytempdir, "/cellranger/"))
        if(!is.null(metricsFilePath)){
          dir.create(paste0(mytempdir, "/cellranger/outs"))
        }
        
        # move files to this sample folder
        file.copy(paste0(dirname(matrixFilePath), "/matrix.mtx.gz"), paste0(mytempdir, "/cellranger/matrix.mtx.gz"))
        file.copy(paste0(dirname(barcodesFilePath), "/barcodes.tsv.gz"), paste0(mytempdir, "/cellranger/barcodes.tsv.gz"))
        file.copy(paste0(dirname(featuresFilePath), "/features.tsv.gz"), paste0(mytempdir, "/cellranger/features.tsv.gz"))
        if(!is.null(metricsFilePath)){
          file.copy(paste0(dirname(metricsFilePath), "/metrics_summary.csv"), paste0(mytempdir, "/cellranger/outs/metrics_summary.csv"))
        }
        
        # make object
        newSce <- importCellRangerV3Sample(dataDir = paste0(mytempdir, "/cellranger/"))
        
        # delete sample folder
        unlink(paste0(mytempdir, "/cellranger/"), recursive = TRUE)
      }
      
    } else if (entry$type == "starSolo") {
      newSce <- importSTARsolo(
        STARsoloDirs = entry$params$STARsoloDirs,
        samples = entry$params$samples,
        STARsoloOuts = entry$params$STARsoloOuts,
        delayedArray = delayedArray
      )
    } else if (entry$type == "starSolo_files") {
      
      # get current tempDir by shiny
      mytempdir <- tempdir(check=TRUE)
      
      if(dir.exists(mytempdir)){
        # get uploaded filepaths
        matrixFilePath <- entry$params$assayFile
        barcodesFilePath <- entry$params$annotFile
        featuresFilePath <- entry$params$featureFile
        # metricsFilePath <- entry$params$summaryFile
        
        # rename to original names
        file.rename(matrixFilePath, paste0(dirname(matrixFilePath), "/matrix.mtx"))
        file.rename(barcodesFilePath, paste0(dirname(barcodesFilePath), "/barcodes.tsv"))
        file.rename(featuresFilePath, paste0(dirname(featuresFilePath), "/features.tsv"))
        # if(!is.null(metricsFilePath)){
        #   file.rename(metricsFilePath, paste0(dirname(metricsFilePath), "/metrics_summary.csv"))
        # }
        
        # create a sample folder
        dir.create(paste0(mytempdir, "/Solo.out/"))
        dir.create(paste0(mytempdir, "/Solo.out/Gene/"))
        dir.create(paste0(mytempdir, "/Solo.out/Gene/filtered"))
        # if(!is.null(metricsFilePath)){
        #   dir.create(paste0(mytempdir, "/cellranger/outs"))
        # }
        
        # move files to this sample folder
        file.copy(paste0(dirname(matrixFilePath), "/matrix.mtx"), paste0(mytempdir, "/Solo.out/Gene/filtered/matrix.mtx"))
        file.copy(paste0(dirname(barcodesFilePath), "/barcodes.tsv"), paste0(mytempdir, "/Solo.out/Gene/filtered/barcodes.tsv"))
        file.copy(paste0(dirname(featuresFilePath), "/features.tsv"), paste0(mytempdir, "/Solo.out/Gene/filtered/features.tsv"))
        # if(!is.null(metricsFilePath)){
        #   file.copy(paste0(dirname(metricsFilePath), "/metrics_summary.csv"), paste0(mytempdir, "/cellranger/outs/metrics_summary.csv"))
        # }
        
        # make object
        newSce <- importSTARsolo(STARsoloDirs = paste0(mytempdir, "/Solo.out"), samples = "sample1")
        
        # delete sample folder
        unlink(paste0(mytempdir, "/Solo.out/"), recursive = TRUE)
        
      # newSce <- importSTARsolo(
      #   STARsoloDirs = entry$params$STARsoloDirs,
      #   samples = entry$params$samples,
      #   STARsoloOuts = entry$params$STARsoloOuts,
      #   delayedArray = delayedArray
      # )
      }
    } else if (entry$type == "busTools") {
      newSce <- importBUStools(
        BUStoolsDirs = entry$params$BUStoolsDirs,
        samples = entry$params$samples,
        delayedArray = delayedArray
      )
    } else if (entry$type == "busTools_files") {
      # get current tempDir by shiny
      mytempdir <- tempdir(check=TRUE)
      
      if(dir.exists(mytempdir)){
        # get uploaded filepaths
        matrixFilePath <- entry$params$assayFile
        barcodesFilePath <- entry$params$annotFile
        featuresFilePath <- entry$params$featureFile
        # metricsFilePath <- entry$params$summaryFile
        
        # rename to original names
        file.rename(matrixFilePath, paste0(dirname(matrixFilePath), "/genes.mtx"))
        file.rename(barcodesFilePath, paste0(dirname(barcodesFilePath), "/genes.barcodes.txt"))
        file.rename(featuresFilePath, paste0(dirname(featuresFilePath), "/genes.genes.txt"))
        # if(!is.null(metricsFilePath)){
        #   file.rename(metricsFilePath, paste0(dirname(metricsFilePath), "/metrics_summary.csv"))
        # }
        
        # create a sample folder
        dir.create(paste0(mytempdir, "/bus_output/"))
        
        # move files to this sample folder
        file.copy(paste0(dirname(matrixFilePath), "/genes.mtx"), paste0(mytempdir, "/bus_output/genes.mtx"))
        file.copy(paste0(dirname(barcodesFilePath), "/genes.barcodes.txt"), paste0(mytempdir, "/bus_output/genes.barcodes.txt"))
        file.copy(paste0(dirname(featuresFilePath), "/genes.genes.txt"), paste0(mytempdir, "/bus_output/genes.genes.txt"))
        
        # make object
        newSce <- importBUStools(BUStoolsDirs = paste0(mytempdir, "/bus_output"), samples = "sample1")
        
        # delete sample folder
        unlink(paste0(mytempdir, "/bus_output/"), recursive = TRUE)
      }
    } else if (entry$type == "seqc") {
      newSce <- importSEQC(
        seqcDirs = entry$params$seqcDirs,
        samples = entry$params$samples,
        prefix = entry$params$prefix,
        delayedArray = delayedArray
      )
    } else if (entry$type == "seqc_files") {
      # get current tempDir by shiny
      mytempdir <- tempdir(check=TRUE)
      
      if(dir.exists(mytempdir)){
        # get uploaded filepaths
        readCountsLocation <- entry$params$readCountsLocation
        moleculeCountsLocation <- entry$params$moleculeCountsLocation
        barcodesLocation <- entry$params$barcodesLocation
        genesLocation <- entry$params$genesLocation
        
        # rename to original names
        file.rename(readCountsLocation, paste0(dirname(readCountsLocation), "/seqc_sparse_read_counts.mtx"))
        file.rename(moleculeCountsLocation, paste0(dirname(moleculeCountsLocation), "/seqc_sparse_molecule_counts.mtx"))
        file.rename(barcodesLocation, paste0(dirname(barcodesLocation), "/seqc_sparse_counts_barcodes.csv"))
        file.rename(genesLocation, paste0(dirname(genesLocation), "/seqc_sparse_counts_genes.csv"))

        # create a sample folder
        dir.create(paste0(mytempdir, "/seqc/"))
        
        # move files to this sample folder
        file.copy(paste0(dirname(readCountsLocation), "/seqc_sparse_read_counts.mtx"), paste0(mytempdir, "/seqc/seqc_sparse_read_counts.mtx"))
        file.copy(paste0(dirname(moleculeCountsLocation), "/seqc_sparse_molecule_counts.mtx"), paste0(mytempdir, "/seqc/seqc_sparse_molecule_counts.mtx"))
        file.copy(paste0(dirname(barcodesLocation), "/seqc_sparse_counts_barcodes.csv"), paste0(mytempdir, "/seqc/seqc_sparse_counts_barcodes.csv"))
        file.copy(paste0(dirname(genesLocation), "/seqc_sparse_counts_genes.csv"), paste0(mytempdir, "/seqc/seqc_sparse_counts_genes.csv"))
        
        # make object
        newSce <- importSEQC(seqcDirs = paste0(mytempdir, "/seqc"), samples = "sample1", prefix = "seqc")
        
        # delete sample folder
        unlink(paste0(mytempdir, "/seqc/"), recursive = TRUE)
      }
    } else if (entry$type == "optimus") {
      newSce <- importOptimus(
        OptimusDirs = entry$params$OptimusDirs,
        samples = entry$params$samples,
        delayedArray = delayedArray
      )
    } else if (entry$type == "optimus_files") {
      # get current tempDir by shiny
      mytempdir <- tempdir(check=TRUE)
      
      if(dir.exists(mytempdir)){
        # get uploaded filepaths
        matrixLocation <- entry$params$matrixLocation
        colIndexLocation <- entry$params$colIndexLocation
        rowIndexLocation <- entry$params$rowIndexLocation
        cellMetricsLocation <- entry$params$cellMetricsLocation
        geneMetricsLocation <- entry$params$geneMetricsLocation
        emptyDropsLocation <- entry$params$emptyDropsLocation
        
        # rename to original names
        file.rename(matrixLocation, paste0(dirname(matrixLocation), "/sparse_counts.npz"))
        file.rename(colIndexLocation, paste0(dirname(colIndexLocation), "/sparse_counts_col_index.npy"))
        file.rename(rowIndexLocation, paste0(dirname(rowIndexLocation), "/sparse_counts_row_index.npy"))
        file.rename(cellMetricsLocation, paste0(dirname(cellMetricsLocation), "/merged-cell-metrics.csv.gz"))
        file.rename(geneMetricsLocation, paste0(dirname(geneMetricsLocation), "/merged-gene-metrics.csv.gz"))
        file.rename(emptyDropsLocation, paste0(dirname(emptyDropsLocation), "/empty_drops_result.csv"))
        
        # create a sample folder
        dir.create(paste0(mytempdir, "/optimus/"))
        dir.create(paste0(mytempdir, "/optimus/call-MergeCountFiles/"))
        dir.create(paste0(mytempdir, "/optimus/call-MergeCellMetrics/"))
        dir.create(paste0(mytempdir, "/optimus/call-MergeGeneMetrics/"))
        dir.create(paste0(mytempdir, "/optimus/call-RunEmptyDrops/"))
        
        # move files to this sample folder
        file.copy(paste0(dirname(matrixLocation), "/sparse_counts.npz"), paste0(mytempdir, "/optimus/call-MergeCountFiles/sparse_counts.npz"))
        file.copy(paste0(dirname(colIndexLocation), "/sparse_counts_col_index.npy"), paste0(mytempdir, "/optimus/call-MergeCountFiles/sparse_counts_col_index.npy"))
        file.copy(paste0(dirname(rowIndexLocation), "/sparse_counts_row_index.npy"), paste0(mytempdir, "/optimus/call-MergeCountFiles/sparse_counts_row_index.npy"))
        file.copy(paste0(dirname(cellMetricsLocation), "/merged-cell-metrics.csv.gz"), paste0(mytempdir, "/optimus/call-MergeCellMetrics/merged-cell-metrics.csv.gz"))
        file.copy(paste0(dirname(geneMetricsLocation), "/merged-gene-metrics.csv.gz"), paste0(mytempdir, "/optimus/call-MergeGeneMetrics/merged-gene-metrics.csv.gz"))
        file.copy(paste0(dirname(emptyDropsLocation), "/empty_drops_result.csv"), paste0(mytempdir, "/optimus/call-RunEmptyDrops/empty_drops_result.csv"))

        # make object
        newSce <- importOptimus(OptimusDirs = paste0(mytempdir, "/optimus"), samples = "sample1")
        
        # delete sample folder
        unlink(paste0(mytempdir, "/optimus/"), recursive = TRUE)
      }
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
    if(entry$type %in% c("cellRanger2", "cellRanger3", "starSolo", "busTools", "seqc", "optimus", "example", "cellRanger3_files", "starSolo_files", "busTools_files", "optimus_files", "seqc_files")){
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

