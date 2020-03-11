#' Apply scMerge batch effect correction method to SingleCellExperiment object
#'
#' The scMerge method leverages factor analysis, stably expressed genes (SEGs) 
#' and (pseudo-) replicates to remove unwanted variations and merge multiple 
#' scRNA-Seq data. 
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param useAssay character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the field 
#' of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"scMerge"`. The name for the corrected 
#' full-sized expression matrix.
#' @param kmeansK vector of int, default `NULL`. A vector indicating the 
#' kmeans' K-value for each batch, in order to construct pseudo-replicates. The 
#' length of `kmeansK` needs to be the same as the number of batches.
#' @param cellType character, default `"cell_type"`. A string indicating the 
#' field of `colData(inSCE)` that defines different cell types.
#' @param seg array, default `NULL`. An array of gene names or indices that 
#' specifies SEG (Stably Expressed Genes) set as negative control. Pre-defined 
#' dataset with human and mouse SEG lists is available to user by running 
#' `data('SEG')`.
#' @param nCores integer, default `parallel::detectCores()`. The number of 
#' cores of processors to allocate for the task. By default it takes all the 
#' cores available to the user. 
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated 
#' with corrected full-sized expression matrix.
#' @export
#' @references Hoa, et al., 2020
#' @examples 
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCMerge(sceBatches, species = 'human')
runSCMerge <- function(inSCE, useAssay = "logcounts", batch = 'batch', 
                       assayName = "scMerge", seg = NULL, kmeansK = NULL, 
                       cellType = 'cell_type', 
                       nCores = parallel::detectCores()){
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
