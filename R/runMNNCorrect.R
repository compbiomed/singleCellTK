#' Apply the mutual nearest neighbors (MNN) batch effect correction method to 
#' SingleCellExperiment object
#' 
#' MNN is designed for batch correction of single-cell RNA-seq data where the 
#' batches are partially confounded with biological conditions of interest. It 
#' does so by identifying pairs of MNN in the high-dimensional log-expression 
#' space. For each MNN pair, a pairwise correction vector is computed by 
#' applying a Gaussian smoothing kernel with bandwidth `sigma`.
#' @param inSCE SingleCellExperiment object. An object that stores your dataset
#' and analysis procedures.
#' @param exprs character, default `"logcounts"`. A string indicating the name 
#' of the assay requiring batch correction in "inSCE", should exist in 
#' `assayNames(inSCE)`.
#' @param batchKey character, default `"batch"`. A string indicating the 
#' field of `colData(inSCE)` that defines different batches.
#' @param reducedDimName character, default `"MNN"`. The name for the 
#' corrected low-dimensional representation.
#' @param nHVG integer, default `1000`. The number of top highly variable genes 
#' to select per batch. Afterwards the intersection of all HVG sets will be 
#' taken for running , thus the exact dimensionality of the resulting
#' reducedDim is not certain. 
#' @param k integer, default `20`. Specifies the number of nearest neighbours to 
#' consider when defining MNN pairs. This should be interpreted as the minimum 
#' frequency of each cell type or state in each batch. Larger values will 
#' improve the precision of the correction by increasing the number of MNN 
#' pairs, at the cost of reducing accuracy by allowing MNN pairs to form between
#' cells of different type.
#' @param sigma Numeric, default `0.1`. Specifies how much information is 
#' shared between MNN pairs when computing the batch effect. Larger values will 
#' share more information, approaching a global correction for all cells in the 
#' same batch. Smaller values allow the correction to vary across cell types, 
#' which may be more accurate but comes at the cost of precision. 
#' @return SingleCellExperiment object with `reducedDim(inSCE, reducedDimName)` 
#' updated with corrected low-dimentional representation.
#' @export
#' @references Lun ATL, et al., 2016 & 2018 
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
#' sceCorr <- runMNNCorrect(sceBatches)
runMNNCorrect <- function(inSCE, exprs = 'logcounts', batchKey = 'batch', 
                          reducedDimName = 'MNN', nHVG = 1000, 
                          k = 20, sigma = 0.1){
    ## Input check
    if(!class(inSCE) == "SingleCellExperiment" && 
       !class(inSCE) == "SCtkExperiment"){
        stop("\"inSCE\" should be a SingleCellExperiment Object.")
    }
    if(!exprs %in% SummarizedExperiment::assayNames(inSCE)) {
        stop(paste("\"exprs\" (assay) name: ", exprs, " not found."))
    }
    if(!batchKey %in% names(SummarizedExperiment::colData(inSCE))){
        stop(paste("\"batchKey name:", batchKey, "not found."))
    }
    reducedDimName <- gsub(' ', '_', reducedDimName)
    
    ## Run algorithm
    # Split the batches
    batches <- list()
    batchCol <- SummarizedExperiment::colData(inSCE)[[batchKey]]
    uniqBatch <- unique(batchCol)
    for(i in uniqBatch){
        batches[[i]] <- inSCE[, batchCol == i]
    }
    
    # Select HVG
    topVarGenesPerBatch <- list()
    for(i in uniqBatch){
        if(nrow(batches[[i]]) <= nHVG){
            topVarGenesPerBatch[[i]] <- 1:nrow(batches[[i]])
        } else {
            mvTrend <- scran::trendVar(batches[[i]], use.spikes=FALSE)
            decomposeTrend <- scran::decomposeVar(batches[[i]], mvTrend)
            topVarGenesPerBatch[[i]] <- order(decomposeTrend$bio, 
                                              decreasing = TRUE)[1:nHVG]
        }    
    }
    selectedHVG <- BiocGenerics::Reduce(intersect, topVarGenesPerBatch)
    
    ## Run mnnCorrect
    inputAssays <- list()
    for(i in uniqBatch){
        inputAssays[[i]] <- 
            SummarizedExperiment::assay(batches[[i]], exprs)[selectedHVG,]
    }
    corrected <- BiocGenerics::do.call(batchelor::mnnCorrect, 
                         c(inputAssays, list(k = k, sigma = sigma)))
    corrected <- t(assay(corrected, 'corrected')[,colnames(inSCE)])
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <- corrected
    
    return(inSCE)
}
