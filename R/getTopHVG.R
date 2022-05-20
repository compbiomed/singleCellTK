#' Get or set top HVG after calculation
#' @description Extracts or select the top variable genes from an input
#' \linkS4class{SingleCellExperiment} object. Note that the variability metrics
#' must be computed using the \code{runFeatureSelection} method before 
#' extracting the feature names of the top variable features. \code{getTopHVG}
#' only returns a character vector of the HVG selection, while with 
#' \code{setTopHVG}, a logical vector of the selection will be saved in the 
#' \code{rowData}, and optionally, a subset object for the HVGs can be stored 
#' in the \code{altExps} slot.
#' @rdname getTopHVG
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param method Specify which method to use for variable gene extraction
#' from Seurat \code{"vst"}, \code{"mean.var.plot"}, \code{"dispersion"} or 
#' Scran \code{"modelGeneVar"}. Default \code{"vst"}
#' @param hvgNumber Specify the number of top variable genes to extract.
#' @param featureDisplay A character string for the \code{rowData} variable name
#' to indicate what type of feature ID should be displayed. If set by 
#' \code{\link{setSCTKDisplayRow}}, will by default use it. If \code{NULL}, will
#' use \code{rownames(inSCE)}.
#' @param altExp \code{TRUE} for also creating a subset \code{inSCE} object with
#' the selected HVGs and store this subset in the \code{altExps} slot, named by
#' \code{rowSubsetName}. Default \code{FALSE}.
#' @param rowSubsetName A character string for the \code{rowData} variable name
#' to store a logical index of selected HVGs. Default 
#' \code{paste0("HVG_", method, hvgNumber)}
#' @return 
#' \item{getTopHVG}{A character vector of the top \code{hvgNumber} variable 
#' feature names}
#' \item{setTopHVG}{The input \code{inSCE} object with the logical vector of 
#' HVG selection updated in \code{rowData}. If \code{altExp} is \code{TRUE},
#' an \code{altExp} is also added}
#' @export
#' @author Irzam Sarfraz, Yichen Wang
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- runSeuratFindHVG(sce)
#' getTopHVG(sce, hvgNumber = 10)
#' sce <- setTopHVG(sce, )
#' @seealso \code{\link{runFeatureSelection}}, \code{\link{runSeuratFindHVG}},
#' \code{\link{runModelGeneVar}}, \code{\link{plotTopHVG}}
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
getTopHVG <- function(inSCE, 
                      method = c("vst", "dispersion", 
                                 "mean.var.plot", "modelGeneVar"), 
                      hvgNumber = 2000, 
                      featureDisplay = metadata(inSCE)$featureDisplay) {
    method <- match.arg(method)
    topGenes <- list()
    if(method == "vst" || method == "dispersion" || method == "modelGeneVar"){
        varianceColumnName = ""
        if (method == "vst") {
            if (is.null(rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized)) {
                stop("'seurat_variableFeatures_vst_varianceStandardized' ",
                     "metric not found in rowData of input sce object. Run ",
                     "Seurat feature selection with 'vst' method before using ",
                     "this function!")
                #inSCE <- runSeuratFindHVG(inSCE, hvgMethod = "vst", 
                #                          hvgNumber = n, altExp = FALSE, ...)
            }
            varianceColumnName = "seurat_variableFeatures_vst_varianceStandardized"
        }
        else if (method == "dispersion") {
            if (is.null(rowData(inSCE)$seurat_variableFeatures_dispersion_dispersion)) {
                stop("'seurat_variableFeatures_dispersion_dispersion' metric ",
                     "not found in rowData of input sce object. Run Seurat ",
                     "feature selection with 'dispersion' method before using ",
                     "this function!")
                #inSCE <- runSeuratFindHVG(inSCE, hvgMethod = "dispersion", 
                #                          hvgNumber = n, altExp = FALSE, ...)
            }
            varianceColumnName = "seurat_variableFeatures_dispersion_dispersion"
        }
        else if (method == "modelGeneVar") {
            if (is.null(rowData(inSCE)$scran_modelGeneVar_bio)) {
                stop("'scran_modelGeneVar_bio' metric not found in rowData of",
                     "input sce object. Run scran feature selection with ",
                     "'modelGeneVar' method before using this function!")
                #inSCE <- scranModelGeneVar(inSCE, ...)
            }
            varianceColumnName = "scran_modelGeneVar_bio"
        }
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE),
            variance = rowData(inSCE)[[varianceColumnName]])
        tempDataFrame <- tempDataFrame[order(-tempDataFrame$variance),]
        tempDataFrame <- tempDataFrame[tempDataFrame$variance > 0, ]
        hvgNumber <- min(hvgNumber, nrow(tempDataFrame))
        topGenes <- as.character(tempDataFrame$featureNames[seq_len(hvgNumber)])
    }
    else if (method == "mean.var.plot") {
        if (is.null(rowData(inSCE)$seurat_variableFeatures_mvp_mean)
            || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersion)
            || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)) {
            stop("'Seurat mean.var.plot' metrics not found in rowData of ",
                 "input sce object. Run Seurat feature selection with ",
                 "'mean.var.plot' method before using this function!")
            #inSCE <- runSeuratFindHVG(inSCE, hvgMethod = "mean.var.plot", 
            #                          hvgNumber = n, altExp = FALSE, ...)
        }
        tempDataFrame <- data.frame(
            featureNames = rownames(inSCE),
            mean = rowData(inSCE)$seurat_variableFeatures_mvp_mean,
            disp = rowData(inSCE)$seurat_variableFeatures_mvp_dispersion,
            dispScaled = rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)
        tempDataFrame <- tempDataFrame[order(-tempDataFrame$disp),]
        tempDataFrame <- tempDataFrame[tempDataFrame$disp > 0,]
        means.use <- (tempDataFrame$mean > 0.1) & (tempDataFrame$mean < 8)
        dispersions.use <- (tempDataFrame$dispScaled > 1) & 
                           (tempDataFrame$dispScaled < Inf)
        tempDataFrame <- tempDataFrame[which(means.use & dispersions.use),]
        hvgNumber <- min(hvgNumber, nrow(tempDataFrame))
        topGenes <- as.character(tempDataFrame$featureNames)[seq_len(hvgNumber)]
    }
    topGenes <- stats::na.omit(topGenes)
    if (!is.null(featureDisplay)) {
        geneIdx <- featureIndex(topGenes, inSCE)
        topGenes <- rowData(inSCE)[[featureDisplay]][geneIdx]
    }
    return(topGenes)
}

#' @rdname getTopHVG
#' @export
#' @importFrom SingleCellExperiment rowSubset
#' @importFrom S4Vectors metadata<-
setTopHVG <- function(inSCE, 
                      method =  c("vst", "dispersion", 
                                  "mean.var.plot", "modelGeneVar"), 
                      hvgNumber = 2000, 
                      rowSubsetName = paste0("HVG_", method, hvgNumber),
                      altExp = FALSE) {
    method <- match.arg(method)
    hvg <- getTopHVG(inSCE, method, hvgNumber, featureDisplay = NULL)
    rowSubset(inSCE, rowSubsetName) <- hvg
    message(paste0(date(), " ... HVG variable '", rowSubsetName, "' created."))
    metadata(inSCE)$sctk$HVG[[rowSubsetName]] <- list(method = method,
                                                      hvgNumber = hvgNumber)
    if (isTRUE(altExp)) {
        inSCE <- subsetSCERows(inSCE, 
                               rowData = rowSubsetName, 
                               returnAsAltExp = TRUE, 
                               altExpName = rowSubsetName)
    }
    return(inSCE)
}
