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
#'  provided \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment.
#' @param inSCE Input \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain a raw counts matrix before empty droplets have been removed.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample separately.
#'  If NULL, then all cells will be processed together. Default NULL.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @param lower See \link[DropletUtils]{emptyDrops} for more information.
#' @param niters See \link[DropletUtils]{emptyDrops} for more information.
#' @param testAmbient See \link[DropletUtils]{emptyDrops} for more information.
#' @param ignore See \link[DropletUtils]{emptyDrops} for more information.
#' @param alpha See \link[DropletUtils]{emptyDrops} for more information.
#' @param retain See \link[DropletUtils]{emptyDrops} for more information.
#' @param barcodeArgs See \link[DropletUtils]{emptyDrops} for more information.
#' @param BPPARAM See \link[DropletUtils]{emptyDrops} for more information.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{emptyDrops} output table appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{emptyDrops_total}, \emph{emptyDrops_logprob},
#'  \emph{emptyDrops_pvalue}, \emph{emptyDrops_limited}, \emph{emptyDrops_fdr}.
#'  Please refer to the documentation of \link[DropletUtils]{emptyDrops} for
#'  details.
#' @examples
#' # The following unfiltered PBMC_1k_v3 data were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' # /pbmc_1k_v3
#' # Only the top 10 cells with most counts and the last 10 cells with non-zero
#' # counts are included in this example.
#' # This example only serves as an proof of concept and a tutorial on how to
#' # run the function. The results should not be
#' # used for drawing scientific conclusions.
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
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))

  
  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }
  
  message(date(), " ... Running 'emptyDrops'")
  
  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
                                 dropletUtils_emptyDrops_total = integer(ncol(inSCE)),
                                 dropletUtils_emptyDrops_logprob = numeric(ncol(inSCE)),
                                 dropletUtils_emptyDrops_pvalue = numeric(ncol(inSCE)),
                                 dropletUtils_emptyDrops_limited = logical(ncol(inSCE)),
                                 dropletUtils_emptyDrops_fdr = numeric(ncol(inSCE)))
  
  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
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
  }
  
  colData(inSCE) = cbind(colData(inSCE), output)
 
  argsList <- argsList[!names(argsList) %in% c("BPPARAM")]
  inSCE@metadata$runEmptyDrops <- argsList[-1]
  inSCE@metadata$runEmptyDrops$packageVersion <- utils::packageDescription("DropletUtils")$Version
  
  return(inSCE)
}
