#' Apply scMerge batch effect correction method to SingleCellExperiment object
#'
#' The scMerge method leverages factor analysis, stably expressed genes (SEGs)
#' and (pseudo-) replicates to remove unwanted variations and merge multiple
#' scRNA-Seq data.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay A single character indicating the name of the assay requiring
#' batch correction. Default \code{"logcounts"}.
#' @param batch A single character indicating a field in
#' \code{\link[SummarizedExperiment]{colData}} that annotates the batches.
#' Default \code{"batch"}.
#' @param kmeansK An integer vector. Indicating the kmeans' K-value for each
#' batch (i.e. how many subclusters in each batch should exist), in order to
#' construct pseudo-replicates. The length of code{kmeansK} needs to be the same
#' as the number of batches. Default \code{NULL}, and this value will be
#' auto-detected by default, depending on \code{cellType}.
#' @param cellType A single character. A string indicating a field in
#' \code{colData(inSCE)} that defines different cell types. Default
#' \code{'cell_type'}.
#' @param seg A vector of gene names or indices that specifies SEG (Stably
#' Expressed Genes) set as negative control. Pre-defined dataset with human and
#' mouse SEG lists is available to user by running \code{data('SEG')}. Default
#' \code{NULL}, and this value will be auto-detected by default with
#' \code{\link[scMerge]{scSEGIndex}}.
#' @param nCores An integer. The number of cores of processors to allocate for
#' the task. Default \code{1L}.
#' @param assayName A single characeter. The name for the corrected assay. Will
#' be saved to \code{\link[SummarizedExperiment]{assay}}. Default
#' \code{"scMerge"}.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{assay(inSCE, assayName)} updated.
#' @export
#' @references Hoa, et al., 2020
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCMerge(sceBatches)
runSCMerge <- function(inSCE, useAssay = "logcounts", batch = 'batch',
                       assayName = "scMerge", seg = NULL, kmeansK = NULL,
                       cellType = 'cell_type',
                       nCores = 1L){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"useAssay\" (assay) name: ", useAssay, " not found."))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    if(is.null(cellType) & is.null(kmeansK)){
        stop("\"cellType\" and \"kmeansK\" cannot be NULL at the same time")
    }
    if(!cellType %in% names(SummarizedExperiment::colData(inSCE))){
        # If NULL, scMerge still works
        stop(paste("\"cellType\" name:", cellType, "not found"))
    }

    nCores <- min(as.integer(nCores), parallel::detectCores())
    assayName <- gsub(' ', '_', assayName)

    ## Run algorithm

    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    uniqBatch <- unique(batchCol)

    # Infer parameters
    if(is.null(cellType)){
        cellTypeCol <- NULL
    } else {
        cellTypeCol <- SummarizedExperiment::colData(inSCE)[[cellType]]
    }
    ## kmeansK
    if(!is.null(cellType) && is.null(kmeansK)){
        # If kmeansK not given, detect by cell type.
        cellTypeCol <- SummarizedExperiment::colData(inSCE)[[cellType]]
        kmeansK <- c()
        for (i in 1:length(uniqBatch)){
            cellTypePerBatch <- cellTypeCol[batchCol == uniqBatch[i]]
            kmeansK <- c(kmeansK, length(unique(cellTypePerBatch)))
        }
        cat("Detected kmeansK:\n")
        print(t(data.frame(K = kmeansK, row.names = uniqBatch)))
    }
    ## SEG
    if(is.null(seg)){
        bpParam <- BiocParallel::MulticoreParam(workers = nCores)
        seg <- scMerge::scSEGIndex(SummarizedExperiment::assay(inSCE, useAssay),
                                   cell_type = cellTypeCol,
                                   BPPARAM = bpParam)
        ctl <- rownames(seg[order(seg$segIdx, decreasing = TRUE)[1:1000],])
    } else {
        ctl <- seg
    }

    # scMerge automatically search for the column called "batch"...
    colDataNames <- names(SummarizedExperiment::colData(inSCE))
    names(SummarizedExperiment::colData(inSCE))[colDataNames == batch] <- 'batch'
    bpParam <- BiocParallel::MulticoreParam(workers = nCores)
    inSCE <- scMerge::scMerge(sce_combine = inSCE, exprs = useAssay,
                              hvg_exprs = useAssay,
                              assay_name = assayName,
                              ctl = ctl, kmeansK = kmeansK,
                              #marker_list = topVarGenesPerBatch,
                              cell_type = cellTypeCol,
                              BPPARAM = bpParam)
    colDataNames <- names(SummarizedExperiment::colData(inSCE))
    names(SummarizedExperiment::colData(inSCE))[colDataNames == 'batch'] <- batch
    return(inSCE)
}
