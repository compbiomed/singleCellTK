#' Apply the mutual nearest neighbors (MNN) batch effect correction method to 
#' SingleCellExperiment object
#'
#' SCANORAMA is analogous to computer vision algorithms for panorama stitching 
#' that identify images with overlapping content and merge these into a larger 
#' panorama. 
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batchKey character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param assayName character, default `"SCANORAMA"`. The name for the 
#' corrected full-sized expression matrix.
#' @param SIGMA numeric, default `15`. Algorithmic parameter, correction 
#' smoothing parameter on Gaussian kernel.
#' @param ALPHA numeric, default `0.1`. Algorithmic parameter, alignment score 
#' minimum cutoff.
#' @param KNN integer, default `20L`. Algorithmic parameter, number of nearest 
#' neighbors to use for matching.
#' @return SingleCellExperiment object with `assay(inSCE, assayName)` updated 
#' with corrected full-sized expression matrix.
#' @export
#' @references Brian Hie et al, 2019
#' @examples 
#' data('sceBatches', package = 'singleCellTK')
#' sceBatches
#' ## class: SingleCellExperiment 
#' ## dim: 27610 1820 
#' ## metadata(0):
#' ## assays(3): normcounts logcounts
#' ## rownames(27610): GCG MALAT1 ... LOC102724004 LOC102724238
#' ## rowData names(0):
#' ## colnames(1820): reads.12732 reads.12733 ... Sample_1598 Sample_1600
#' ## colData names(2): cell_type1 batch
#' ## reducedDimNames(5): PCA
#' ## spikeNames(0):
#' sceCorr <- runSCANORAMA(sceBatches)
runSCANORAMA <- function(inSCE, exprs = 'logcounts', batchKey = 'batch', 
                         assayName = 'SCANORAMA', SIGMA = 15, ALPHA = 0.1, 
                         KNN = 20L){
    ## Input check
    if(!class(inSCE) == "SingleCellExperiment" && 
       !class(inSCE) == "SCtkExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found"))
    }
    if(!batchKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batchKey\" name:", batchKey, "not found"))
    }
    assayName <- gsub(' ', '_', assayName)
    
    ## Run algorithm
    batchCol <- SummarizedExperiment::colData(inSCE)[[batchKey]]
    batches <- unique(batchCol)
    nBatch <- length(batches)
    ixList <- list()
    for(i in 1:nBatch){
        ixList[[i]] <- which(batchCol == batches[[i]])
    }
    exprsList <- list()
    for(i in 1:nBatch){
        exprsList[[i]] <- t(as.matrix(SummarizedExperiment::assay(inSCE[,ixList[[i]]], exprs)))
    }
    geneList <- list()
    for(i in 1:nBatch){
        geneList[[i]] <- colnames(exprsList[[i]])
    }
    
    #TODO: add reducedDim support
    corr <- scnrm$correct(datasets_full = exprsList, 
                          genes_list = geneList,   
                          sigma = SIGMA, alpha = ALPHA, 
                          knn = KNN, verbose = 0)
    
    ## corr[[1]] for corrected expressions
    ##    corr[[1]][[n]] corr-exprs of the n-th batch, in Python 
    ##    scipy.sparse.csr.scr_matrix class
    ## corr[[2]] The gene list
    
    corrExpList <- list()
    for(i in 1:nBatch){
        corrExpList[[i]] <- corr[[1]][[i]]$toarray()
    }
    
    #TODO: Find better language for this concatenation
    outMatrix <- rbind(corrExpList[[1]], corrExpList[[2]])
    if(nBatch > 2){
        for(i in 3:nBatch){
            outMatrix <- rbind(outMatrix, corrExpList[[i]])
        }
    } 
    
    # Note that in the scnrm output the ordering is "shuffled"
    # Then we'll have to order it back to match the original SCE
    
    #TODO: Find better language for this ordering
    geneOrder <- c()
    for(i in 1:nrow(inSCE)){
        geneOrder <- c(geneOrder, which(rownames(inSCE)[i] == corr[[2]]))
    }
    #TODO: Find better language for this concatenation
    cellOrder <- c(ixList[[1]], ixList[[2]])
    if(nBatch > 2){
        for(i in 3:nBatch){
            cellOrder <- c(cellOrder, ixList[[i]])
        }
    }
    
    outMatrix <- t(outMatrix[cellOrder, geneOrder])
    SummarizedExperiment::assay(inSCE, assayName) <- outMatrix
    return(inSCE)
}
