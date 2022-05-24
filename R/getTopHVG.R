#' Get or set top HVG after calculation
#' @description Extracts or select the top variable genes from an input
#' \linkS4class{SingleCellExperiment} object. Note that the variability metrics
#' must be computed using the \code{runFeatureSelection} method before 
#' extracting the feature names of the top variable features. \code{getTopHVG}
#' only returns a character vector of the HVG selection, while with 
#' \code{setTopHVG}, a logical vector of the selection will be saved in the 
#' \code{rowData}, and optionally, a subset object for the HVGs can be stored 
#' in the \code{altExps} slot at the same time.
#' @rdname getTopHVG
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param method Specify which method to use for variable gene extraction
#' from Seurat \code{"vst"}, \code{"mean.var.plot"}, \code{"dispersion"} or 
#' Scran \code{"modelGeneVar"}. Default \code{"vst"}
#' @param hvgNumber Specify the number of top variable genes to extract.
#' @param useHVGList For getting top HVG, use the HVG list set by 
#' \code{setTopHVG}. Will ignore \code{method} and \code{hvgNumber} if not 
#' \code{NULL}. Default \code{NULL}.
#' @param featureDisplay A character string for the \code{rowData} variable name
#' to indicate what type of feature ID should be displayed. If set by 
#' \code{\link{setSCTKDisplayRow}}, will by default use it. If \code{NULL}, will
#' use \code{rownames(inSCE)}.
#' @param genes A customized character vector of gene list to be set as a 
#' \code{rowData} variable. Will ignore \code{method} and \code{hvgNumber} if 
#' set. Default \code{NULL}.
#' @param genesBy If setting customized \code{genes}, where should it be found
#' in \code{rowData}? Leave \code{NULL} for matching \code{rownames}. Default 
#' \code{NULL}.
#' @param altExp \code{TRUE} for also creating a subset \code{inSCE} object with
#' the selected HVGs and store this subset in the \code{altExps} slot, named by
#' \code{hvgListName}. Default \code{FALSE}.
#' @param hvgListName A character string for the \code{rowData} variable name
#' to store a logical index of selected HVGs. Default \code{NULL}, will be 
#' determined basing on other parameters. 
#' @return 
#' \item{getTopHVG}{A character vector of the top \code{hvgNumber} variable 
#' feature names}
#' \item{setTopHVG}{The input \code{inSCE} object with the logical vector of 
#' HVG selection updated in \code{rowData}, and related parameter updated in 
#' \code{metadata}. If \code{altExp} is \code{TRUE}, an \code{altExp} is also 
#' added}
#' @export
#' @author Irzam Sarfraz, Yichen Wang
#' @examples
#' data("scExample", package = "singleCellTK")
#' sce <- runSeuratFindHVG(sce)
#' hvgs <- getTopHVG(sce, hvgNumber = 10)
#' sce <- setTopHVG(sce, method = "vst", hvgNumber = 5)
#' @seealso \code{\link{runFeatureSelection}}, \code{\link{runSeuratFindHVG}},
#' \code{\link{runModelGeneVar}}, \code{\link{plotTopHVG}}
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
getTopHVG <- function(inSCE, 
                      method = c("vst", "dispersion", 
                                 "mean.var.plot", "modelGeneVar"), 
                      hvgNumber = 2000, 
                      useHVGList = NULL,
                      featureDisplay = metadata(inSCE)$featureDisplay) {
    method <- match.arg(method)
    topGenes <- character()
    if (!is.null(useHVGList)) {
        topGenes <- rownames(inSCE)[rowSubset(inSCE, useHVGList)]
    } else {
        if (method == "vst" || method == "dispersion" || method == "modelGeneVar") {
            varianceColumnName = ""
            if (method == "vst") {
                if (is.null(rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized)) {
                    stop("'seurat_variableFeatures_vst_varianceStandardized' ",
                         "metric not found in rowData of input sce object. ",
                         "Run `runSeuratFindHVG()` with 'vst' method before ", 
                         "using this function!")
                }
                varianceColumnName = "seurat_variableFeatures_vst_varianceStandardized"
            } else if (method == "dispersion") {
                if (is.null(rowData(inSCE)$seurat_variableFeatures_dispersion_dispersion)) {
                    stop("'seurat_variableFeatures_dispersion_dispersion' ", 
                         "metric not found in rowData of input sce object. ", 
                         "Run `runSeuratFindHVG()` with 'dispersion' method ", 
                         "before using this function!")
                }
                varianceColumnName = "seurat_variableFeatures_dispersion_dispersion"
            } else if (method == "modelGeneVar") {
                if (is.null(rowData(inSCE)$scran_modelGeneVar_bio)) {
                    stop("'scran_modelGeneVar_bio' metric not found in ", 
                         "rowData of input sce object. Run ",
                         "`runModelGeneVar()` before using this function!")
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
        } else if (method == "mean.var.plot") {
            if (is.null(rowData(inSCE)$seurat_variableFeatures_mvp_mean)
                || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersion)
                || is.null(rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled)) {
                stop("'Seurat mean.var.plot' metrics not found in rowData of ",
                     "input sce object. Run `runSeuratFindHVG()` with ",
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
                      genes = NULL, genesBy = NULL,
                      hvgListName = NULL,
                      altExp = FALSE) {
    method <- match.arg(method)
    if (!is.null(genes)) {
        # Set customized gene list
        if (is.null(genesBy)) {
            hvg <- genes
        } else {
            geneIdx <- featureIndex(features = genes, inSCE = inSCE, 
                                    by = genesBy)
            hvg <- rownames(inSCE)[geneIdx]
        }
    } else {
        # Use pre-calculated variability metrics
        hvg <- getTopHVG(inSCE, method = method, hvgNumber = hvgNumber, 
                         featureDisplay = NULL)
    }
    if (is.null(hvgListName)) {
        if (!is.null(genes)) {
            hvgListName <- "CustomizedHVG"
        } else {
            hvgListName <- paste0("HVG_", method, hvgNumber)
        }
    }
    rowSubset(inSCE, hvgListName) <- hvg
    message(paste0(date(), " ... HVG variable '", hvgListName, "' created."))
    useAssay <- metadata(inSCE)$sctk$runFeatureSelection[[method]]$useAssay
    metadata(inSCE)$sctk$hvgLists[[hvgListName]] <- list(method = method,
                                                         hvgNumber = hvgNumber,
                                                         useAssay = useAssay)
    if (isTRUE(altExp)) {
        inSCE <- subsetSCERows(inSCE, 
                               rowData = hvgListName, 
                               returnAsAltExp = TRUE, 
                               altExpName = hvgListName)
    }
    return(inSCE)
}
