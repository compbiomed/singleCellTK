#' getTopHVG
#' Extracts the top variable genes from an input \code{SingleCellExperiment} object.
#' Note that the variability metrics must be computed using the `runFeatureSelection`
#' method before extracting the feature names of the top variable features. If
#' `altExp` parameter is a \code{character} value, this function will return the
#' input \code{SingleCellExperiment} object with the subset containing only the
#' top variable features stored as an \code{altExp} slot in returned object.
#' However, if this parameter is set to \code{NULL}, only the names of the top
#' variable features will be returned as a \code{character} vector.
#' @param inSCE Input \code{SingleCellExperiment} object
#' @param method Specify which method to use for variable gene extraction
#' from either Seurat "vst", "mean.var.plot", "dispersion" or Scran
#' "modelGeneVar".
#' @param n Specify the number of top variable genes to extract.
#' @param altExp A \code{character} value that specifies the name of the \code{altExp}
#'  slot that should be created to store the subset \code{SingleCellExperiment}
#'  object containing only the top `n` variable features. Default value is
#'  \code{NULL}, which will not store the subset \code{SingleCellExperiment} object
#'  and instead will only return the names of the top `n` variable features.
#' @return A \code{character} vector of the top variable feature names or the
#'  input \code{SingleCellExperiment} object with subset of variable features
#'  stored as an \code{altExp} in the object. 
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scranModelGeneVar(sce_chcl, "counts")
#' # return top 10 variable genes
#' topGenes <- getTopHVG(sce_chcl, "modelGeneVar", 10)
#' @importFrom SummarizedExperiment rowData
getTopHVG <- function(inSCE, method, n = 2000, altExp = NULL) {
    topGenes <- list()
    if(method == "vst" || method == "dispersion" || method == "modelGeneVar"){
        varianceColumnName = ""
        if (method == "vst") {
            if (is.null(rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized)) {
                stop("'seurat_variableFeatures_vst_varianceStandardized' ",
                     "metric not found in rowData of input sce object. Run ",
                     "Seurat feature selection with 'vst' method before using ",
                     "this function!")
            }
            varianceColumnName = "seurat_variableFeatures_vst_varianceStandardized"
        }
        else if (method == "dispersion") {
            if (is.null(rowData(inSCE)$seurat_variableFeatures_dispersion_dispersion)) {
                stop("'seurat_variableFeatures_dispersion_dispersion' metric ",
                     "not found in rowData of input sce object. Run Seurat ",
                     "feature selection with 'dispersion' method before using ",
                     "this function!")
            }
            varianceColumnName = "seurat_variableFeatures_dispersion_dispersion"
        }
        else if (method == "modelGeneVar") {
            if (is.null(rowData(inSCE)$scran_modelGeneVar_bio)) {
                stop("'scran_modelGeneVar_bio' metric not found in rowData of",
                     "input sce object. Run scran feature selection with ",
                     "'modelGeneVar' method before using this function!")
            }
            varianceColumnName = "scran_modelGeneVar_bio"
        }
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE),
            variance = rowData(inSCE)[varianceColumnName])
        tempDataFrame <-
          tempDataFrame[order(-tempDataFrame[varianceColumnName]),]
        
        tempDataFrame <- 
            tempDataFrame[tempDataFrame[varianceColumnName] > 0, ]
        
        if(nrow(tempDataFrame) < n){
            n <- nrow(tempDataFrame)
        }
        
        topGenes <- as.character(tempDataFrame$featureNames[seq_len(n)])
        
        topGenes <- stats::na.omit(topGenes)
        
        if(!is.null(altExp)){
          topGenes <- inSCE[topGenes, ]
          expData(inSCE = inSCE, assayName = altExp, tag = "hvg", altExp = TRUE) <- topGenes
          topGenes <- inSCE
        }
        
    }
    else if (method == "mean.var.plot") {
        if (is.null(rowData(inSCE)$seurat_variableFeatures_mvp_mean)
            || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersion)
            || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)) {
            stop("'Seurat mean.var.plot' metrics not found in rowData of ",
                 "input sce object. Run Seurat feature selection with ",
                 "'mean.var.plot' method before using this function!")
        }
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE),
            mean = rowData(inSCE)$seurat_variableFeatures_mvp_mean,
            disp = rowData(inSCE)$seurat_variableFeatures_mvp_dispersion,
            dispScaled = rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)
        tempDataFrame <- tempDataFrame[order(-tempDataFrame$disp),]
        
        tempDataFrame <- 
            tempDataFrame[tempDataFrame["disp"] > 0, ]
        
        if(nrow(tempDataFrame) < n){
            n <- nrow(tempDataFrame)
        }
        
        means.use <- (tempDataFrame[, "mean"] > 0.1) &
                     (tempDataFrame[, "mean"] < 8)
        dispersions.use <- (tempDataFrame[, "dispScaled"] > 1) &
                           (tempDataFrame[, "dispScaled"] < Inf)
        topGenes <- as.character(tempDataFrame$featureNames[which(x = means.use & dispersions.use)])[seq_len(n)]
        
        topGenes <- stats::na.omit(topGenes)
        
        if(!is.null(altExp)){
            topGenes <- inSCE[topGenes, ]
            expData(inSCE = inSCE, assayName = altExp, tag = "hvg", altExp = TRUE) <- topGenes
            topGenes <- inSCE
        }
    }
    
    
    return(topGenes)
}
