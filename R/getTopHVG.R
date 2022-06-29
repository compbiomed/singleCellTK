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
#' @param useFeatureSubset Get the feature names in the HVG list set by
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
#' @param featureSubsetName A character string for the \code{rowData} variable
#' name to store a logical index of selected features. Default \code{NULL}, will
#' be determined basing on other parameters.
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
                      useFeatureSubset = NULL,
                      featureDisplay = metadata(inSCE)$featureDisplay) {
    method <- match.arg(method)
    topGenes <- character()
    if (!is.null(useFeatureSubset)) {
        topGenes <- .parseUseFeatureSubset(inSCE, useFeatureSubset,
                                           altExpObj = NULL, returnType = "cha")
    } else {
        metrics <- .dfFromHVGMetric(inSCE, method)
        metrics <- metrics[order(-metrics$v_rank),]
        metrics <- metrics[metrics$v_rank > 0, ]
        if (method == "mean.var.plot") {
            means.use <- (metrics$mean > 0.1) & (metrics$mean < 8)
            dispersions.use <- (metrics$v_plot > 1) & (metrics$v_plot < Inf)
            metrics <- metrics[means.use & dispersions.use,]
        }
        hvgNumber <- min(hvgNumber, nrow(metrics))
        topGenes <- as.character(metrics$featureNames)[seq_len(hvgNumber)]
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
                      featureSubsetName = NULL,
                      genes = NULL, genesBy = NULL,
                      altExp = FALSE) {
    method <- match.arg(method)
    features <- character()
    useAssay <- NULL
    if (!is.null(genes)) {
        # Set customized gene list
        features <- genes
        if (!is.null(genesBy)) {
            geneIdx <- featureIndex(features = features, inSCE = inSCE,
                                    by = genesBy)
            features <- rownames(inSCE)[geneIdx]
        }
    } else {
        # Use pre-calculated variability metrics
        features <- getTopHVG(inSCE, method = method, hvgNumber = hvgNumber,
                              featureDisplay = NULL)
        useAssay <- metadata(inSCE)$sctk$runFeatureSelection[[method]]$useAssay
    }
    nFeature <- length(features)
    if (is.null(featureSubsetName)) {
        if (!is.null(genes)) {
            featureSubsetName <- "CustomizedHVG"
        } else {
            featureSubsetName <- paste0("HVG_", method, hvgNumber)
        }
    }
    rowSubset(inSCE, featureSubsetName) <- features
    message(paste0(date(), " ... Feature subset variable '", featureSubsetName,
                   "' created."))
    metadata(inSCE)$sctk$featureSubsets[[featureSubsetName]] <-
        list(method = method,
             number = nFeature,
             useAssay = useAssay)
    if (isTRUE(altExp)) {
        inSCE <- subsetSCERows(inSCE,
                               rowData = featureSubsetName,
                               returnAsAltExp = TRUE,
                               altExpName = featureSubsetName)
    }
    return(inSCE)
}

#' Check if specified method has already been performed, and extract the metric.
#' Will be used by `getTopHVG` to rank the top, and `plotTopHVG` for the axis
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param method Specify which method to use for variable gene extraction
#' from Seurat \code{"vst"}, \code{"mean.var.plot"}, \code{"dispersion"} or
#' Scran \code{"modelGeneVar"}. Default \code{"vst"}
#' @return data.frame object of HVG metrics calculated by \code{method},
#' containing columns of \code{"mean"}, \code{"v_rank"}, \code{"v_plot"}
#' @noRd
.dfFromHVGMetric <- function(inSCE,
                             method = c("vst", "mean.var.plot", "dispersion",
                                        "modelGeneVar")) {
    method <- match.arg(method)
    df <- data.frame(featureNames = rownames(inSCE))
    if (method == "vst") {
        m <- "seurat_variableFeatures_vst_mean"
        v_rank <- "seurat_variableFeatures_vst_varianceStandardized"
        v_plot <- "seurat_variableFeatures_vst_varianceStandardized"
        if (is.null(rowData(inSCE)[[v_rank]]) || is.null(rowData(inSCE)[[m]])) {
            stop("Seurat vst metric not found in inSCE. ",
                 "Run `runSeuratFindHVG()` with 'vst' method before ",
                 "using this function!")
        }
    } else if (method == "dispersion") {
        m <- "seurat_variableFeatures_dispersion_mean"
        v_rank <- "seurat_variableFeatures_dispersion_dispersion"
        v_plot <- "seurat_variableFeatures_dispersion_dispersionScaled"
        if (is.null(rowData(inSCE)[[v_rank]]) ||
            is.null(rowData(inSCE)[[m]]) ||
            is.null(rowData(inSCE)[[v_plot]])) {
            stop("Seurat dispersion metric not found in inSCE. ",
                 "Run `runSeuratFindHVG()` with 'dispersion' method ",
                 "before using this function!")
        }
    } else if (method == "modelGeneVar") {
        m <- "scran_modelGeneVar_mean"
        v_rank <- "scran_modelGeneVar_bio"
        v_plot <- "scran_modelGeneVar_totalVariance"
        if (is.null(rowData(inSCE)[[v_rank]]) ||
            is.null(rowData(inSCE)[[m]]) ||
            is.null(rowData(inSCE)[[v_plot]])) {
            stop("Scran modelGeneVar metric not found in inSCE. Run ",
                 "`runModelGeneVar()` before using this function!")
        }
    } else if (method == "mean.var.plot") {
        m <- "seurat_variableFeatures_mvp_mean"
        v_rank <- "seurat_variableFeatures_mvp_dispersion"
        v_plot <- "seurat_variableFeatures_mvp_dispersionScaled"
        if (is.null(rowData(inSCE)[[v_rank]]) ||
            is.null(rowData(inSCE)[[m]]) ||
            is.null(rowData(inSCE)[[v_plot]])) {
            stop("Seurat mean.var.plot metric not found in inSCE. ",
                 "Run `runSeuratFindHVG()` with ",
                 "'mean.var.plot' method before using this function!")
        }
    }
    df$mean <- rowData(inSCE)[[m]]
    df$v_rank <- rowData(inSCE)[[v_rank]]
    df$v_plot <- rowData(inSCE)[[v_plot]]
    return(df)
}

#' parse `useFeatureSubset` in other functions such as `scaterPCA`, `getUMAP`..
#' Do checks, and return logical vector. Or character vector as needed by Seurat
#' methods
#' @param inSCE Input \linkS4class{SingleCellExperiment} object
#' @param useFeatureSubset Subset of feature to use. A character string
#' indicating a \code{rowData} variable that stores the logical vector of HVG
#' selection, or a vector that can subset the rows of \code{inSCE}.
#' @param altExpObj A \linkS4class{SingleCellExperiment} object, extracted from
#' \code{altExps} of \code{inSCE}, can be \code{identical()} to \code{inSCE}.
#' Used when functions like \code{\link{scaterPCA}} allows \code{useAltExp}, so
#' that the output vector match to the rownames of \code{altExpObj}
#' @param returnType What type of vector should be returned
#' @return \code{NULL} if \code{is.null(useFeatureSubset)}. Otherwise, a vector
#' subsetting the features.
#' @noRd
.parseUseFeatureSubset <- function(inSCE, useFeatureSubset, altExpObj = NULL,
                                  returnType = c("logical", "character")) {
    returnType <- match.arg(returnType)
    if (!is.null(altExpObj)) {
        if (!inherits(altExpObj, "SingleCellExperiment")) {
            stop()
        }
    }
    if (identical(altExpObj, inSCE)) altExpObj <- NULL
    features <- NULL
    if (!is.null(useFeatureSubset)) {
        if (is.character(useFeatureSubset) & length(useFeatureSubset) == 1) {
            # Assume using rowData variable
            if (!useFeatureSubset %in% names(rowData(inSCE))) {
                stop("Specified feature subset not found.")
            }
            feature.logical <- rowData(inSCE)[[useFeatureSubset]]
            if (!is.logical(feature.logical)) {
                stop("Specified rowData variable is not logical.")
            }
            features <- rownames(inSCE)[feature.logical]
            if (useFeatureSubset %in% names(metadata(inSCE)$sctk$featureSubsets)) {
                # If the subset is created with SCTK HVG methods
                # Return a ranked gene list
                method <- metadata(inSCE)$sctk$featureSubsets[[useFeatureSubset]]$method
                useAssay <- metadata(inSCE)$sctk$featureSubsets[[useFeatureSubset]]$useAssay
                metricAssay <- metadata(inSCE)$sctk$runFeatureSelection[[method]]$useAssay
                if (useAssay != metricAssay) {
                    warning("Sorting features subset using metrics ",
                            "calculated from an assay different than the ",
                            "assay used for the subset creation")
                }
                metrics <- .dfFromHVGMetric(inSCE, method)
                metrics <- metrics[order(-metrics$v_rank),]
                features <- metrics$featureNames[metrics$featureNames %in% features]
            }
        } else if (is.character(useFeatureSubset) &
                   length(useFeatureSubset) > 1) {
            # Assume using rownames
            if (any(!useFeatureSubset %in% rownames(inSCE))) {
                nf <- useFeatureSubset[!useFeatureSubset %in% rownames(inSCE)]
                stop("Specified features not found: ",
                     paste(nf, collapse = ", "))
            }
            features <- useFeatureSubset
        } else if (is.logical(useFeatureSubset)) {
            if (length(useFeatureSubset) != nrow(inSCE)) {
                stop("Length of logical subset doesn't match nrow(inSCE).")
            }
            features <- rownames(inSCE)[useFeatureSubset]
        } else if (is.numeric(useFeatureSubset)) {
            features <- rownames(inSCE)[useFeatureSubset]
            if (any(is.na(features))) {
                stop("Numeric subscript contains out-of-bounds indices.")
            }
        } else {
            stop("Invalid useFeatureSubset specification.")
        }
        if (!is.null(altExpObj)) {
            features <- features[features %in% rownames(altExpObj)]
        }
        if (returnType == "logical") {
            if (!is.null(altExpObj)) {
                features <- rownames(altExpObj) %in% features
            } else {
                features <- rownames(inSCE) %in% features
            }
        }

    }
    return(features)
}
