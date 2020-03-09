#' Apply scMerge batch effect correction method to SingleCellExperiment object
#'
#' The scMerge method leverages factor analysis, stably expressed genes (SEGs) 
#' and (pseudo-) replicates to remove unwanted variations and merge multiple 
#' scRNA-Seq data. 
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batch character, default `"batch"`. A string indicating the field 
#' of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"scMerge"`. The name for the corrected 
#' full-sized expression matrix.
#' @param kmeansK vector of int, default `NULL`. Indicates the kmeans's K for 
#' each batch. The length of `kmeansK` needs to be the same as the number of 
#' batches. If not given, this vector will be identified by counting cell types 
#' in each batch.
#' @param cellType character, default `"cell_type"`. A string indicating the 
#' field of `colData(inSCE)` that defines different cell types. Only needed 
#' when `kmeansK` is left to `NULL`. 
#' @param species character, default `NULL`. Choose from `{"human", "mouse"}`. 
#' If given, the algorithm will take default species specific SEG (Stably 
#' Expressed Genes) set as negative control. 
#' @param seg character vecter, default `NULL`. A gene list that specifies SEG 
#' (Stably Expressed Genes) set as negative control. If `species` set, will be 
#' ignored; if both `species` and `seg` are `NULL`, will take some time to 
#' automatically detect top 1000 SEGs from `assay(inSCE, exprs)`.
#' @param nHVG integer, default `1000`. The number of top highly variable genes 
#' to select per batch.
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated 
#' with corrected full-sized expression matrix.
#' @export
#' @references Hoa, et al., 2020
#' @examples 
#' data('sceBatches', package = 'singleCellTK')
#' sceCorr <- runSCMerge(sceBatches, species = 'human')
runSCMerge <- function(inSCE, exprs = "logcounts", batch = 'batch', 
                       assayName = "scMerge", kmeansK = NULL, 
                       cellType = 'cell_type', species = NULL, 
                       seg = NULL, nHVG = 1000){
    ## Input check
    if(!inherits(inSCE, "SingleCellExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found."))
    }
    if(!batch %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batch\" name:", batch, "not found"))
    }
    if(!cellType %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"cellType\" name:", cellType, "not found"))
    }
    if (!is.null(species) && !is.null(seg)){
        stop("None or only one of the arguments \"species\" and \"seg\" should be applied.")
    } 
    
    assayName <- gsub(' ', '_', assayName)
    
    if(!is.null(species)){
        if(species %in% c('human', 'mouse')){
            # Automatic selection
            genes <- rownames(inSCE)
            if(all(startsWith(genes, 'ENG'))){
                message("All gene IDs start with \"ENG\", using Ensembl IDs")
                ctl <- SEG[[species]][[1]]
            } else {
                message("Not all gene IDs start with \"ENG\", using gene symbols")
                ctl <- SEG[[species]][[2]]
            }
        } else {
            stop("Please choose \"species\" from {\"human\", \"mouse\"}")
        }
    } else if (!is.null(seg)){
        ctl <- seg
    } else {
        seg <- scMerge::scSEGIndex(SummarizedExperiment::assay(inSCE, exprs), 
                            SummarizedExperiment::colData(inSCE)[[cellType]],
                            ncore = parallel::detectCores())
        ctl <- rownames(seg[order(seg$segIdx, decreasing = TRUE)[1:1000],])
    }
    
    ## Run algorithm
    
    # Select HVG
    batches <- list()
    batchCol <- SummarizedExperiment::colData(inSCE)[[batch]]
    uniqBatch <- unique(batchCol)
    for(i in uniqBatch){
        batches[[i]] <- inSCE[, batchCol == i]
    }
    topVarGenesPerBatch <- list()
    for(i in uniqBatch){
        if(nrow(batches[[i]]) <= nHVG){
            topVarGenesPerBatch[[i]] <- 1:nrow(batches[[i]])
        } else {
          sce.var <- scran::modelGeneVar(batches[[i]])
          topVarGenesPerBatch[[i]] <- order(sce.var$bio,decreasing = TRUE)[seq(nHVG)]
            #mvTrend <- scran::trendVar(batches[[i]], use.spikes=FALSE)
            #decomposeTrend <- scran::decomposeVar(batches[[i]], mvTrend)
            #topVarGenesPerBatch[[i]] <- order(decomposeTrend$bio, 
                                              #decreasing = TRUE)[1:nHVG]
        }    
    }
    
    if(is.null(kmeansK)){
        # If kmeansK not given, detect by cell type.
        cellTypeCol <- SummarizedExperiment::colData(inSCE)[[cellType]]
        kmeansK <- c()
        for (i in 1:length(uniqBatch)){
            cellTypePerBatch <- cellTypeCol[batchCol == uniqBatch[i]]
            kmeansK <- c(kmeansK, length(unique(cellTypePerBatch)))
        }
        print("Detected kmeansK:")
        print(kmeansK)
    }
    
    # scMerge automatically search for the column called "batch"...
    sceTmp <- inSCE
    colDataNames <- names(SummarizedExperiment::colData(sceTmp))
    names(SummarizedExperiment::colData(sceTmp))[colDataNames == batch] <- 'batch'
    sceTmp <- scMerge::scMerge(sceTmp, exprs = exprs, 
                               hvg_exprs = exprs, 
                               assay_name = assayName, 
                               ctl = ctl, kmeansK = kmeansK, 
                               marker_list = topVarGenesPerBatch)
    corrected <- SummarizedExperiment::assay(sceTmp, assayName)
    SummarizedExperiment::assay(inSCE, assayName) <- corrected
    return(inSCE)
}
