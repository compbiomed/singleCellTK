#' @title Detecting contamination with DecontX.
#' @description A wrapper function for \link[celda]{decontX}. Identify
#'  potential contamination from experimental factors such as ambient RNA.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample A single character specifying a name that can be found in
#' \code{colData(inSCE)} to directly use the cell annotation; or a character
#' vector with as many elements as cells to indicates which sample each cell
#' belongs to. Default NULL. \link[celda]{decontX} will be run on cells from
#' each sample separately.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#' 'counts'.
#' @param background A \link[SingleCellExperiment]{SingleCellExperiment}
#' with the matrix located in the assay slot under \code{bgAssayName}. It should have 
#' the same structure as inSCE except it contains the matrix of empty droplets instead 
#' of cells. When supplied, empirical distribution of transcripts from these 
#' empty droplets will be used as the contamination distribution. Default NULL.
#' @param bgAssayName Character. Name of the assay to use if background is a 
#' \link[SingleCellExperiment]{SingleCellExperiment}. If NULL, the function
#' will use the same value as \code{useAssay}. Default is NULL. 
#' @param bgBatch Batch labels for \code{background}. If \code{background} is a 
#' \link[SingleCellExperiment]{SingleCellExperiment} object, this can be a single 
#' character specifying a name that can be found in \code{colData(background)} 
#' to directly use the barcode annotation; or a numeric / character vector that has  
#' as many elements as barcodes to indicate which sample each barcode belongs to. Its 
#' unique values should be the same as those in \code{sample}, such that each 
#' batch of cells have their corresponding batch of empty droplets as background, 
#' pointed by this parameter. Default to NULL.
#' @param z Numeric or character vector. Cell cluster labels. If NULL,
#' PCA will be used to reduce the dimensionality of the dataset initially,
#' '\link[uwot]{umap}' from the 'uwot' package
#' will be used to further reduce the dataset to 2 dimenions and
#' the '\link[dbscan]{dbscan}' function from the 'dbscan' package
#' will be used to identify clusters of broad cell types. Default NULL.
#' @param maxIter Integer. Maximum iterations of the EM algorithm. Default 500.
#' @param convergence Numeric. The EM algorithm will be stopped if the maximum
#' difference in the contamination estimates between the previous and
#' current iterations is less than this. Default 0.001.
#' @param iterLogLik Integer. Calculate log likelihood every \code{iterLogLik}
#' iteration. Default 10.
#' @param delta Numeric Vector of length 2. Concentration parameters for
#' the Dirichlet prior for the contamination in each cell. The first element
#' is the prior for the native counts while the second element is the prior for
#' the contamination counts. These essentially act as pseudocounts for the
#' native and contamination in each cell. If \code{estimateDelta = TRUE},
#' this is only used to produce a random sample of proportions for an initial
#' value of contamination in each cell. Then
#' \code{\link[MCMCprecision]{fit_dirichlet}} is used to update
#' \code{delta} in each iteration.
#' If \code{estimateDelta = FALSE}, then \code{delta} is fixed with these
#' values for the entire inference procedure. Fixing \code{delta} and
#' setting a high number in the second element will force \code{decontX}
#' to be more aggressive and estimate higher levels of contamination at
#' the expense of potentially removing native expression.
#' Default \code{c(10, 10)}.
#' @param estimateDelta Boolean. Whether to update \code{delta} at each
#' iteration.
#' @param varGenes Integer. The number of variable genes to use in
#' dimensionality reduction before clustering. Variability is calcualted using
#' \code{\link[scran]{modelGeneVar}} function from the 'scran' package.
#' Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter
#' used in '\link[dbscan]{dbscan}' to estimate broad cell clusters.
#' Used only when z is not provided. Default 1.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  'decontX_Contamination' and 'decontX_Clusters' added to the
#'  \link{colData} slot. Additionally, the
#' decontaminated counts will be added as an assay called 'decontXCounts'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDecontX(sce[,sample(ncol(sce),20)])
#' @export
runDecontX <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    background = NULL,
    bgAssayName = NULL,
    bgBatch = NULL,
    z = NULL,
    maxIter = 500,
    delta = c(10, 10),
    estimateDelta = TRUE,
    convergence = 0.001,
    iterLogLik = 10,
    varGenes = 5000,
    dbscanEps = 1,
    seed = 12345,
    logfile = NULL,
    verbose = TRUE
) {
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  if(!is.null(sample)) {
    # get sample from colData(inSCE)
    if (length(sample) == 1) {
      if (!sample %in% names(SummarizedExperiment::colData(inSCE))) {
        stop("Specified Sample variable not found in colData")
      }
      sample <- SummarizedExperiment::colData(inSCE)[[sample]]
    } else if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
    if (is.factor(sample)) {
      sample <- as.character(sample)
    }
    uniqueSample <- unique(sample)

    if (!is.null(background)) {
      ### Background must be a SCE object if the input of count is in SCE object. Required by celda::decontX. 
      if (!inherits(background, 'SingleCellExperiment')) {
        stop("'background' is not a SingleCellExperiment object")
      }

      if (!is.null(bgBatch)) {
        ### check whether variable can be found in colData.
        if(length(bgBatch) == 1) { # & is(background, 'SingleCellExperiment')
          if(!bgBatch %in% names(SummarizedExperiment::colData(background))) {
            stop("Specified bgBatch variable not found in colData")
          }
          bgBatch <- SummarizedExperiment::colData(background)[[bgBatch]]
        } else if(length(bgBatch) != ncol(background)) {
          stop("'bgBatch' must be the same length as the number of columns in 'background'")      
        }

        if (!all(unique(bgBatch) %in% uniqueSample)) {
            stop("Not all bgBatch can be found in 'samples'.")
        }

        if (is.factor(bgBatch)) {
          bgBatch <- as.character(bgBatch)
        }
      }      
    }
  } else {
    # When sample not specified, labeled as "all_cells" 
    sample <- rep("all_cells", ncol(inSCE))
    if (!is.null(background)) {
        if (!is.null(bgBatch)) {
            warning("Using all background because 'sample' not specified.")
        }
        bgBatch <- rep("all_cells", ncol(background))
    }
    uniqSample <- "all_cells"    
  }



  message(paste0(date(), " ... Running 'DecontX'"))

  ## keep a copy of original decontX result, if there's any
  oldColData <- NULL
  oldColId <- grep('decontX_', colnames(SummarizedExperiment::colData(inSCE)))
  if(length(oldColId) != 0) {
    oldColData <- SummarizedExperiment::colData(inSCE)[, oldColId]
    SummarizedExperiment::colData(inSCE)[oldColId] <- NULL
  }

  oldAssay <- NULL
  oldAssayId <- grep('decontXcounts', SummarizedExperiment::assayNames(inSCE))
  if(length(oldAssayId) != 0) {
    oldAssay <- SummarizedExperiment::assays(inSCE)[oldAssayId]
    SummarizedExperiment::assays(inSCE)[oldAssayId] <- NULL
  }

  oldDim<- NULL
  oldDimId <- grep('decontX_', SingleCellExperiment::reducedDimNames(inSCE))
  if(length(oldAssayId) != 0) {
    oldDim <- SingleCellExperiment::reducedDims(inSCE)[oldDimId]
    SingleCellExperiment::reducedDims(inSCE)[oldDimId] <- NULL
  }

  oldMeta <- NULL
  oldMetaId <- grep('runDecontX', names(S4Vectors::metadata(inSCE)$sctk))
  if(length(oldMetaId) != 0) {
    oldMeta <- S4Vectors::metadata(inSCE)$sctk[oldMetaId]
    S4Vectors::metadata(inSCE)$sctk[oldMetaId] <- NULL
  }

  rm.ix <- which(colSums(assay(inSCE, useAssay)) == 0)
  if(length(rm.ix) > 0){
    inSCEOrig <- inSCE
    inSCE <- inSCE[,-rm.ix]
    sample <- sample[-rm.ix]
  }

  inSCE <- celda::decontX(x = inSCE,
                          batch = sample,
                          assayName = useAssay,
                          background = background,
                          bgAssayName = bgAssayName,
                          bgBatch = bgBatch,
                          z = z,
                          maxIter = maxIter,
                          delta = delta,
                          estimateDelta = estimateDelta,
                          convergence = convergence,
                          iterLogLik = iterLogLik,
                          varGenes = varGenes,
                          dbscanEps = dbscanEps,
                          seed = seed,
                          logfile = logfile,
                          verbose = verbose)

  #argsList <- argsList[!names(argsList) %in% ("...")]
  colId <- colnames(SummarizedExperiment::colData(inSCE)) %in% c('decontX_contamination', 'decontX_clusters')
  newCol <- SummarizedExperiment::colData(inSCE)[colId]
  SummarizedExperiment::colData(inSCE)[colId] <- NULL

  assayId <- SummarizedExperiment::assayNames(inSCE) == "decontXcounts"
  newAssay <- SummarizedExperiment::assays(inSCE)[assayId]
  SummarizedExperiment::assays(inSCE)[assayId] <- NULL

  DimId <- grep("decontX", SingleCellExperiment::reducedDimNames(inSCE), 
                value=TRUE)#SingleCellExperiment::reducedDimNames(inSCE) == "decontX_UMAP"
  newDim <- SingleCellExperiment::reducedDims(inSCE)[DimId]
  SingleCellExperiment::reducedDims(inSCE)[DimId] <- NULL

  metaId <- names(S4Vectors::metadata(inSCE)) == "decontX"
  newMeta <- S4Vectors::metadata(inSCE)[metaId]
  newMeta$decontX <- S4Vectors::metadata(inSCE)$decontX$runParams #only keep runParams in metadata
  newMeta$decontX$packageVersion <- utils::packageDescription("celda")$Version
  names(newMeta) <- "runDecontX"
  S4Vectors::metadata(inSCE)[metaId] <- NULL

  ### further process the metadata for each sample
  newMeta_bc <- newMeta
  newMeta <- NULL
  batchMeta <- newMeta_bc$runDecontX$batch
  bgBatchMeta <- newMeta_bc$runDecontX$batchBackground
  newMeta_bc$runDecontX[c("batch", "batchBackground")] <- NULL

  for (s in unique(sample)) {
    newMeta$runDecontX[[s]] <- newMeta_bc$runDecontX
    newMeta$runDecontX[[s]]$batch <- batchMeta[batchMeta == s]
    newMeta$runDecontX[[s]]$batchBackground <- bgBatchMeta[bgBatchMeta == s]
  }

  ### merge inSCE to restore some cells that are removed by rm.ix
  if(length(rm.ix) > 0){
    inSCE <- mergeSCEColData(inSCE1 = inSCEOrig, inSCE2 = inSCE)
  }

  ### update the decontX result, if there's any
  
  if(!is.null(background)) {
    ## need to rename the colData if backgroud is supplied
    colnames(newCol) <- paste(colnames(newCol), 'bg', sep = '_')
    ## need to rename the assay name
    names(newAssay) <- "decontXcounts_bg"
    ## need to rename reducedDim
    names(newDim) <- paste(names(newDim), "bg", sep = "_")
    ## need to rename metadata
    names(newMeta) <- "runDecontX_bg"
  }

  # update Assay
  if(!is.null(oldAssay)) {
    oldAssay[names(newAssay)] <- newAssay
    finalAssay <- oldAssay
  } else {
    finalAssay <- newAssay
  }

  # update reducedDim
  if(!is.null(oldDim)) {
    oldDim[names(newDim)] <- newDim
    finalDim <- oldDim
  } else {
    finalDim <- newDim
  }

  # update colData
  if(!is.null(oldColData)) {
    oldColData[, colnames(newCol)] <- newCol[rownames(oldColData), ]
    finalCol <- oldColData    
  } else {
    finalCol <- newCol[colnames(inSCE), ]
  }

  # update metadata
  if(!is.null(oldMeta)) {
    oldMeta[names(newMeta)] <- newMeta
    finalMeta <- oldMeta
  } else {
    finalMeta <- newMeta
  }
  
  # store new result back to inSCE object
  SummarizedExperiment::colData(inSCE)[colnames(finalCol)] <- finalCol
  SummarizedExperiment::assays(inSCE)[names(finalAssay)] <- finalAssay
  SingleCellExperiment::reducedDims(inSCE)[names(finalDim)] <- finalDim
  S4Vectors::metadata(inSCE)$sctk[names(finalMeta)] <- finalMeta
  # inSCE@metadata$sctk$runDecontX <- argsList[-1]
  # inSCE@metadata$sctk$runDecontX$packageVersion <- utils::packageDescription("celda")$Version
  inSCE <- expSetDataTag(inSCE, "raw", "decontXcounts")
  return(inSCE)
}

