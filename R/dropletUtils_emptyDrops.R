.runEmptyDrops <- function(barcode.matrix=barcode.matrix, lower=lower,
                           niters=niters,
                           test.ambient=test.ambient,
                           ignore=ignore, 
                           alpha=alpha,
                           retain=retain,
                           barcode.args=list(),
                           BPPARAM=BiocParallel::SerialParam()) {
  
  barcode.matrix <- .convertToMatrix(barcode.matrix)
  
  result <- DropletUtils::emptyDrops(m = barcode.matrix, 
                                     lower = lower,
                                     niters = niters,
                                     test.ambient = test.ambient,
                                     ignore = ignore,
                                     alpha = alpha,
                                     retain = retain,
                                     barcode.args = barcode.args,
                                     BPPARAM = BPPARAM)
  colnames(result) <- paste0("dropletUtils_emptyDrops_", colnames(result))
  
  return(result)
}


#' @title Identify empty droplets using \link[DropletUtils]{emptyDrops}.
#' @description Run \link[DropletUtils]{emptyDrops} on the count matrix in the
#' provided \linkS4class{SingleCellExperiment} object.
#' Distinguish between droplets containing cells and ambient RNA in a
#' droplet-based single-cell RNA sequencing experiment.
#' @param inSCE A \linkS4class{SingleCellExperiment} object. Must contain a raw 
#' counts matrix before empty droplets have been removed.
#' @param sample Character vector or colData variable name. Indicates which 
#' sample each cell belongs to. Default \code{NULL}.
#' @param useAssay A string specifying which assay in the SCE to use. Default 
#' \code{"counts"}
#' @param lower See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{100}.
#' @param niters See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{10000}.
#' @param testAmbient See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{FALSE}. 
#' @param ignore See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{NULL}.
#' @param alpha See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{NULL}.
#' @param retain See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{NULL}.
#' @param barcodeArgs See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{list()}.
#' @param BPPARAM See \link[DropletUtils]{emptyDrops} for more information. 
#' Default \code{BiocParallel::SerialParam()}.
#' @return A \linkS4class{SingleCellExperiment} object with the
#' \link[DropletUtils]{emptyDrops} output table appended to the
#' \link{colData} slot. The columns include
#' \code{emptyDrops_total}, \code{emptyDrops_logprob},
#' \code{emptyDrops_pvalue}, \code{emptyDrops_limited}, \code{emptyDrops_fdr}.
#' Please refer to the documentation of \link[DropletUtils]{emptyDrops} for
#' details.
#' @seealso \code{\link{runDropletQC}}, \code{\link{plotEmptyDropsResults}},
#' \code{\link{plotEmptyDropsScatter}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- runEmptyDrops(inSCE = sce)
#' @import DropletUtils
#' @export
#' @importFrom SummarizedExperiment colData colData<-
runEmptyDrops <- function(inSCE,
                          sample = NULL,
                          useAssay = "counts", 
                          lower = 100,
                          niters = 10000,
                          testAmbient = FALSE,
                          ignore = NULL, 
                          alpha = NULL,
                          retain = NULL,
                          barcodeArgs = list(),
                          BPPARAM = BiocParallel::SerialParam()
) {
  # getting the current argument values
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  argsList <- argsList[!names(argsList) %in% c("inSCE", "BPPARAM")]
  argsList$packageVersion <- utils::packageDescription("DropletUtils")$Version
  
  sample <- .manageCellVar(inSCE, var = sample)
  if (is.null(sample)) {
    sample <- rep(1, ncol(inSCE))
  }
  
  message(date(), " ... Running 'emptyDrops'")
  
  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(
    row.names = colnames(inSCE),
    dropletUtils_emptyDrops_total = integer(ncol(inSCE)),
    dropletUtils_emptyDrops_logprob = numeric(ncol(inSCE)),
    dropletUtils_emptyDrops_pvalue = numeric(ncol(inSCE)),
    dropletUtils_emptyDrops_limited = logical(ncol(inSCE)),
    dropletUtils_emptyDrops_fdr = numeric(ncol(inSCE))
  )
  
  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (s in samples) {
    sceSampleInd <- sample == s
    sceSample <- inSCE[, sceSampleInd]
    
    mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
    result <- .runEmptyDrops(barcode.matrix = mat,
                             lower = lower,
                             niters = niters,
                             test.ambient = testAmbient,
                             ignore = ignore, 
                             alpha = alpha,
                             retain = retain,
                             barcode.args = barcodeArgs,
                             BPPARAM = BPPARAM)
    
    
    output[sceSampleInd, ] <- result
    S4Vectors::metadata(output[sceSampleInd, ]) <- S4Vectors::metadata(result)
    if (!identical(samples, 1)) {
      S4Vectors::metadata(inSCE)$sctk$runEmptyDrops[[s]] <- argsList
    }
  }
  if (identical(samples, 1)) {
    S4Vectors::metadata(inSCE)$sctk$runEmptyDrops$all_cells <- argsList
  }
  
  colData(inSCE) <- cbind(colData(inSCE), output)
  return(inSCE)
}
