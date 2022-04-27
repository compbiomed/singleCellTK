.matrixTypeCheck <- function(inSCE, redDimType = NULL,
                             useAssay = NULL, useReducedDim = NULL,
                             useAltExp = NULL) {
  # Helper function for checking if the specified matrix type is valid
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  if (redDimType == "embedding") {
    if (is.null(useAssay) && is.null(useReducedDim)) {
      stop("`useAssay` and `useReducedDim` cannot be NULL at the same time.")
    } else if (!is.null(useAssay) && !is.null(useReducedDim)) {
      stop("`useAssay` and `useReducedDim` cannot be specified at the same time.")
    } else {
      if (!is.null(useReducedDim)) {
        if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
          stop("Specified `useReducedDim` not found.")
        }
        if (!is.null(useAltExp)) {
          warning("`useAltExp` will be ignored when using `useReducedDim`.")
        }
        sce <- inSCE
      } else {
        if (!is.null(useAltExp)) {
          if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
            stop("Specified `useAltExp` not found.")
          }
          sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
          if (!useAssay %in% SummarizedExperiment::assayNames(sce)) {
            stop("Specified `useAssay` not found in `useAltExp`.")
          }
        } else {
          if (!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
            stop("Specified `useAssay` not found.")
          }
          sce <- inSCE
        }
      }
    }
  } else if (redDimType == "linear") {
    if (!is.null(useReducedDim)) {
      stop("Currently `useReducedDim` is not allowed for linear dimension reduction.")
    }
    if (is.null(useAssay)) {
      stop("`useAssay` cannot be NULL for linear dimension reduction")
    } else {
      if (is.null(useAltExp)) {
        if (!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
          stop("Specified `useAssay` not found.")
        }
      } else {
        if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
          stop("Specified `useAltExp` not found.")
        }
        sce <- SingleCellExperiment::altExp(inSCE, useAltExp)
        if (!useAssay %in% SummarizedExperiment::assayNames(sce)) {
          stop("Specified `useAssay` not found in `useAltExp`.")
        }
      }
    }
  }
}

#' Generic Wrapper function for running dimensionality reduction
#' @details Wrapper function to run one of the available dimensionality
#' reduction algorithms integrated within SCTK from \code{\link{scaterPCA}},
#' \code{\link{runSeuratPCA}}, \code{\link{runSeuratICA}}, \code{\link{getTSNE}},
#' \code{\link{runSeuratTSNE}}, \code{\link{getUMAP}} and
#' \code{\link{runSeuratUMAP}}. Users can use an assay by specifying
#' \code{useAssay}, use the assay in an altExp by specifying both
#' \code{useAltExp} and \code{useAssay}, or use a low-dimensionality
#' representation by specifying \code{useReducedDim}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param method One from \code{"scaterPCA"}, \code{"seuratPCA"},
#' \code{"seuratICA"}, \code{"rTSNE"}, \code{"seuratTSNE"}, \code{"scaterUMAP"}
#' and \code{"seuratUMAP"}.
#' @param useAssay Assay to use for computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"counts"}.
#' @param useAltExp The subset to use for computation, usually for the
#' selected variable features. Default \code{NULL}.
#' @param useReducedDim The low dimension representation to use for embedding
#' computation. Default \code{NULL}.
#' @param reducedDimName The name of the result matrix. Required.
#' @param nComponents Specify the number of dimensions to compute with the
#'  selected method in case of PCA/ICA and the number of components to
#'  use in the case of TSNE/UMAP methods.
#' @param seed Random seed for reproducibility of results.
#' Default \code{NULL} will use global seed in use by the R environment.
#' @param ... The other arguments for running a specific algorithm. Please refer
#' to the one you use.
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{reducedDim} updated with the result.
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runNormalization(sce, useAssay = "counts",
#'                         outAssayName = "logcounts_scaled",
#'                         normalizationMethod = "logNormCounts",
#'                         scale = TRUE)
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA",
#'                     useAssay = "logcounts_scaled",
#'                     reducedDimName = "PCA")
runDimReduce <- function(inSCE,
                         method = c("scaterPCA",
                                    "seuratPCA",
                                    "seuratICA",
                                    "rTSNE",
                                    "seuratTSNE",
                                    "scaterUMAP",
                                    "seuratUMAP"),
                         useAssay = NULL, useReducedDim = NULL,
                         useAltExp = NULL, reducedDimName, nComponents = 20, seed = NULL, ...
) {

  method <- match.arg(method)
  args <- list(...)
  if (is.null(reducedDimName)) {
    stop("Must specify `reducedDimName` to store the result.")
  }
  if (method %in% c("scaterPCA", "seuratPCA", "seuratICA")) {
    .matrixTypeCheck(inSCE, "linear", useAssay, useReducedDim, useAltExp)
  } else {
    .matrixTypeCheck(inSCE, "embedding", useAssay, useReducedDim, useAltExp)
  }

  if (method == "scaterPCA") {
    message(paste0(date(), " ... Computing Scater PCA."))
    inSCE <- scaterPCA(inSCE = inSCE, useAssay = useAssay, useAltExp = useAltExp,
                       reducedDimName = reducedDimName, nComponents = nComponents, seed = seed, ...)
  } else if (method == "scaterUMAP") {
    message(paste0(date(), " ... Computing Scater UMAP."))
    inSCE <- getUMAP(inSCE = inSCE, useAssay = useAssay, useAltExp = useAltExp,
                     useReducedDim = useReducedDim,
                     reducedDimName = reducedDimName, seed = seed, ...)
  } else if (method == "rTSNE") {
    message(paste0(date(), " ... Computing RtSNE."))
    inSCE <- getTSNE(inSCE = inSCE, useAssay = useAssay, useAltExp = useAltExp,
                     useReducedDim = useReducedDim,
                     reducedDimName = reducedDimName, seed = seed, ...)
  } else {
    # Seurat part
    if (!is.null(useAltExp)) {
      tempSCE <- SingleCellExperiment::altExp(inSCE, useAltExp)
      # tempSCE <- runSeuratFindHVG(inSCE = tempSCE, useAssay = useAssay,
      #                          altExp = TRUE)
    } else if (!is.null(useAssay)) {
      tempSCE <- inSCE
      #tempSCE <- runSeuratFindHVG(inSCE = tempSCE, useAssay = useAssay)
    }
    if (method %in% c("seuratPCA", "seuratICA")) {
      ## SeuratPCA/ICA
      if (method == "seuratPCA") {
        message(paste0(date(), " ... Computing Seurat PCA."))
        tempSCE <- runSeuratPCA(tempSCE, useAssay = useAssay,
                             reducedDimName = reducedDimName,
                             nPCs = nComponents, features = rownames(inSCE), seed = seed, ...)
      } else if (method == "seuratICA") {
        message(paste0(date(), " ... Computing Seurat ICA."))
        tempSCE <- runSeuratICA(tempSCE, useAssay = useAssay,
                             reducedDimName = reducedDimName,
                             nics = nComponents, seed = seed, ...)
      }
      seuratObj <- tempSCE@metadata$seurat
      if (!is.null(useAltExp)) {
        altExp(inSCE, useAltExp)@metadata$seurat <- seuratObj
      } else if (!is.null(useAssay)) {
        inSCE@metadata$seurat <- seuratObj
      }
    } else {
      ## SeuratUMAP/TSNE
      if (is.null(useReducedDim)) {
        ### using assay
        if (!"useReduction" %in% names(args)) {
          stop("Must specify `useReduction` when using `useAssay` in seuratUMAP/TSNE")
        }
        if (args$useReduction == "pca") {
          message(paste0(date(), " ... Computing Seurat PCA."))
          tempSCE <- runSeuratPCA(inSCE = tempSCE,
                               useAssay = useAssay,
                               reducedDimName = paste0(useAssay, "_seuratPCA"), seed = seed)
        } else if (args$useReduction == "ica") {
          message(paste0(date(), " ... Computing Seurat ICA."))
          tempSCE <- runSeuratICA(inSCE = tempSCE,
                               useAssay = useAssay,
                               reducedDimName = paste0(useAssay, "_seuratICA"), seed = seed)
        }
        if (method == "seuratUMAP") {
          message(paste0(date(), " ... Computing Seurat UMAP."))
          tempSCE <- runSeuratUMAP(inSCE = tempSCE,
                                   reducedDimName = reducedDimName, seed = seed, ...)
        } else {
          message(paste0(date(), " ... Computing Seurat tSNE."))
          tempSCE <- runSeuratTSNE(inSCE = tempSCE,
                                   reducedDimName = reducedDimName, seed = seed, ...)
        }
      } else {
        ### using external reducedDim
        if (!is.null(args$useReduction)) {
          stop("Cannot specify `useReduction` when using `useReducedDim` in seuratUMAP/TSNE")
        }
        tempSCE <- inSCE
        seuratObj <- convertSCEToSeurat(inSCE)
        tempSCE@metadata$seurat$obj <- seuratObj
        reDim <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)
        colnames(reDim) <- paste0(useReducedDim, "_", seq_len(length(colnames(reDim))))
        rownames(reDim) <- gsub('_', '-', rownames(reDim))
        key <-  gsub('_', '', useReducedDim)
        # hard-code "pca"
        tempSCE@metadata$seurat$obj@reductions$pca <-
          Seurat::CreateDimReducObject(embeddings = reDim,
                                       key = paste0(key, "_"), assay = "RNA")
        if (method == "seuratUMAP") {
          # hard-code useReduction="pca"
          message(paste0(date(), " ... Computing Seurat UMAP."))
          tempSCE <- runSeuratUMAP(inSCE = tempSCE, useReduction = "pca",
                                   reducedDimName = reducedDimName, seed = seed, ...)
        } else {
          # hard-code useReduction="pca"
          message(paste0(date(), " ... Computing Seurat tSNE."))
          tempSCE <- runSeuratTSNE(inSCE = tempSCE, useReduction = "pca",
                                   reducedDimName = reducedDimName, seed = seed, ...)
        }
      }
    }
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <-
      SingleCellExperiment::reducedDim(tempSCE, reducedDimName)
  }
  return(inSCE)
}
