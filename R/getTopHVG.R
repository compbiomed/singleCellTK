#' getTopHVG
#' Extracts the top variable genes from an input singleCellExperiment object
#' @param inSCE an input singleCellExperiment object
#' @param method represents which method to use for variable gene extraction from
#' either Seurat "vst", "mean.var.plot", "dispersion" or Scran "modelGeneVar"
#' @param n number of top variable genes to extract
#' @return list of top variable gene names
#' @export
#' @author Irzam Sarfraz
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- scran_modelGeneVar(sce_chcl, "counts")
#' topGenes <- getTopHVG(sce_chcl, "modelGeneVar", 10) #return top 10 variable genes
#' print(topGenes)

getTopHVG <- function(inSCE, method, n = 2000) {
    topGenes <- list()
    if(method == "vst" || method == "dispersion" || method == "modelGeneVar"){
        varianceColumnName = ""
        if(method == "vst"){
            varianceColumnName = "seurat_variableFeatures_vst_varianceStandardized"
        }
        else if(method == "dispersion"){
            varianceColumnName = "seurat_variableFeatures_dispersion_dispersion"
        }
        else if(method == "modelGeneVar"){
            varianceColumnName = "scran_modelGeneVar_bio"
        }
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE), 
            variance = rowData(inSCE)[varianceColumnName])
        tempDataFrame <- tempDataFrame[order(-tempDataFrame[varianceColumnName]),]
        topGenes <- as.character(tempDataFrame$featureNames[1:n])
    }
    else if(method == "mean.var.plot"){
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE), 
            mean = rowData(inSCE)$seurat_variableFeatures_mvp_mean, 
            disp = rowData(inSCE)$seurat_variableFeatures_mvp_dispersion, 
            dispScaled = rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)
        tempDataFrame <- tempDataFrame[order(-tempDataFrame$disp),]
        means.use <- (tempDataFrame[, "mean"] > 0.1) & (tempDataFrame[, "mean"] < 8)
        dispersions.use <- (tempDataFrame[, "dispScaled"] > 1) & (tempDataFrame[, "dispScaled"] < Inf)
        topGenes <- as.character(tempDataFrame$featureNames[which(x = means.use & dispersions.use)])[1:n]
    }
    return(topGenes)
}