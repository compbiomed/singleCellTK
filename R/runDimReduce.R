#' Generic Wrapper function for running dimensionality reduction
#' @details Wrapper function to run one of the available dimensionality
#' reduction algorithms integrated within SCTK from \code{\link{scaterPCA}},
#' \code{\link{runSeuratPCA}}, \code{\link{runSeuratICA}}, \code{\link{runTSNE}},
#' \code{\link{runSeuratTSNE}}, \code{\link{runUMAP}} and
#' \code{\link{runSeuratUMAP}}. Users can use an assay by specifying
#' \code{useAssay}, use the assay in an altExp by specifying both
#' \code{useAltExp} and \code{useAssay}, or use a low-dimensionality
#' representation by specifying \code{useReducedDim}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param method One from \code{"scaterPCA"}, \code{"seuratPCA"},
#' \code{"seuratICA"}, \code{"rTSNE"}, \code{"seuratTSNE"}, \code{"scaterUMAP"},
#' \code{"seuratUMAP"}, \code{"scanpyPCA"}, \code{"scanpyUMAP"} and \code{"scanpyTSNE"}.
#' @param useAssay Assay to use for computation. If \code{useAltExp} is
#' specified, \code{useAssay} has to exist in
#' \code{assays(altExp(inSCE, useAltExp))}. Default \code{"counts"}.
#' @param useAltExp The subset to use for computation, usually for the
#' selected variable features. Default \code{NULL}.
#' @param useReducedDim The low dimension representation to use for embedding
#' computation. Default \code{NULL}.
#' @param reducedDimName The name of the result matrix. Required.
#' @param useFeatureSubset Subset of feature to use for dimension reduction. A
#' character string indicating a \code{rowData} variable that stores the logical
#' vector of HVG selection, or a vector that can subset the rows of
#' \code{inSCE}. Default \code{NULL}.
#' @param scale Logical scalar, whether to standardize the expression values.
#' Default \code{TRUE}.
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
#'                         outAssayName = "logcounts",
#'                         normalizationMethod = "logNormCounts")
#' sce <- runDimReduce(inSCE = sce, method = "scaterPCA",
#'                     useAssay = "logcounts", scale = TRUE,
#'                     reducedDimName = "PCA")
runDimReduce <- function(inSCE,
                         method = c("scaterPCA",
                                    "seuratPCA",
                                    "seuratICA",
                                    "scanpyPCA",
                                    "rTSNE",
                                    "seuratTSNE",
                                    "scaterUMAP",
                                    "seuratUMAP",
                                    "scanpyUMAP",
                                    "scanpyTSNE"),
                         useAssay = NULL, useReducedDim = NULL,
                         useAltExp = NULL, reducedDimName = method,
                         nComponents = 20, useFeatureSubset = NULL,
                         scale = FALSE, seed = 12345, ...)
{

  method <- match.arg(method)
  args <- list(...)
  if (method %in% c("scaterPCA", "seuratPCA", "seuratICA") &
      !is.null(useReducedDim)) {
    stop("`useReducedDim` is not allowed for linear dimension reduction.")
  }

  if (method == "scaterPCA") {
    inSCE <- scaterPCA(inSCE = inSCE, useAssay = useAssay,
                       useAltExp = useAltExp, reducedDimName = reducedDimName,
                       nComponents = nComponents,
                       useFeatureSubset = useFeatureSubset, scale = scale,
                       seed = seed, ...)
  } else if (method == "scaterUMAP") {
    inSCE <- runUMAP(inSCE = inSCE, useAssay = useAssay, useAltExp = useAltExp,
                     useReducedDim = useReducedDim, initialDims = 25,
                     useFeatureSubset = useFeatureSubset, scale = scale,
                     reducedDimName = reducedDimName, seed = seed, ...)
  } else if (method == "scanpyPCA"){
    inSCE <- runScanpyPCA(inSCE = inSCE,
                          useAssay = useAssay, 
                          reducedDimName = reducedDimName, 
                          nPCs = nComponents, 
                          method = "auto", 
                          use_highly_variable = FALSE
                          )
  } else if (method == "scanpyTSNE"){
    inSCE <- runScanpyTSNE(inSCE = inSCE, useAssay = useAssay,
                           useReducedDim = useReducedDim, reducedDimName = reducedDimName, ...)
  } else if (method == "scanpyUMAP"){
    inSCE <- runScanpyUMAP(inSCE = inSCE, useAssay = useAssay, 
                           useReducedDim = useReducedDim, reducedDimName = reducedDimName, ...)
  } else if (method == "rTSNE") {
    inSCE <- runTSNE(inSCE = inSCE, useAssay = useAssay, useAltExp = useAltExp,
                     useReducedDim = useReducedDim,
                     useFeatureSubset = useFeatureSubset, scale = scale,
                     reducedDimName = reducedDimName, seed = seed, ...)
  } else {
    # Seurat part
    # TODO: Honestly, the input checks should have been implemented for
    # functions being wrapped because they are being exposed to users as well.
    # We should not being performing redundant checks when wrapping them again.
    useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                               useReducedDim = useReducedDim,
                               useAltExp = useAltExp, returnMatrix = FALSE)
    useAssay <- useMat$names$useAssay
    if (!is.null(useAltExp)) {
      tempSCE <- SingleCellExperiment::altExp(inSCE, useAltExp)
    } else if (!is.null(useAssay)) {
      tempSCE <- inSCE
    }
    if (method %in% c("seuratPCA", "seuratICA")) {
      ## SeuratPCA/ICA
      if (method == "seuratPCA") {
        p <- paste0(date(), " ... Computing Seurat PCA.")
        message(p)
        tempSCE <- runSeuratPCA(tempSCE, useAssay = useAssay,
                                reducedDimName = reducedDimName,
                                nPCs = nComponents,
                                useFeatureSubset = useFeatureSubset,
                                scale = scale, seed = seed, ...)
      } else if (method == "seuratICA") {
        p <- paste0(date(), " ... Computing Seurat ICA.")
        message(p)
        tempSCE <- runSeuratICA(tempSCE, useAssay = useAssay,
                                reducedDimName = reducedDimName,
                                nics = nComponents,
                                useFeatureSubset = useFeatureSubset,
                                scale = scale, seed = seed, ...)
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
          p <- paste0(date(), " ... Computing Seurat PCA.")
          message(p)
          tempSCE <- runSeuratPCA(inSCE = tempSCE,
                               useAssay = useAssay,
                               reducedDimName = paste0(useAssay, "_seuratPCA"),
                               useFeatureSubset = useFeatureSubset, seed = seed)
        } else if (args$useReduction == "ica") {
          p <- paste0(date(), " ... Computing Seurat ICA.")
          message(p)
          tempSCE <- runSeuratICA(inSCE = tempSCE,
                               useAssay = useAssay,
                               reducedDimName = paste0(useAssay, "_seuratICA"),
                               useFeatureSubset = useFeatureSubset, seed = seed)
        }
        if (method == "seuratUMAP") {
          p <- paste0(date(), " ... Computing Seurat UMAP.")
          message(p)
          tempSCE <- runSeuratUMAP(inSCE = tempSCE,
                                   reducedDimName = reducedDimName,
                                   seed = seed, ...)
        } else {
          p <- paste0(date(), " ... Computing Seurat tSNE.")
          message(p)
          tempSCE <- runSeuratTSNE(inSCE = tempSCE,
                                   reducedDimName = reducedDimName,
                                   seed = seed, ...)
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
          p <- paste0(date(), " ... Computing Seurat UMAP.")
          message(p)
          tempSCE <- runSeuratUMAP(inSCE = tempSCE, useReduction = "pca",
                                   reducedDimName = reducedDimName,
                                   seed = seed, ...)
        } else {
          # hard-code useReduction="pca"
          p <- paste0(date(), " ... Computing Seurat tSNE.")
          message(p)
          tempSCE <- runSeuratTSNE(inSCE = tempSCE, useReduction = "pca",
                                   reducedDimName = reducedDimName,
                                   seed = seed, ...)
        }
      }
    }
    SingleCellExperiment::reducedDim(inSCE, reducedDimName) <-
      SingleCellExperiment::reducedDim(tempSCE, reducedDimName)
  }
  return(inSCE)
}
