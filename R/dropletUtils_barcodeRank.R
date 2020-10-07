
.runBarcodeRankDrops <- function(barcode.matrix, lower=lower,
                                 fit.bounds=fit.bounds,
                                 df=df) {
  
  ## Convert to sparse matrix if not already in that format
  barcode.matrix <- .convertToMatrix(barcode.matrix)
  
  output <- DropletUtils::barcodeRanks(m = barcode.matrix, lower=lower,
                                       fit.bounds=fit.bounds,
                                       df=df)
  
  knee.ix <- as.integer(output@listData$total >= S4Vectors::metadata(output)$knee)
  inflection.ix <- as.integer(output@listData$total >= S4Vectors::metadata(output)$inflection)
  rank.ix<- as.integer(output$rank)
  total.ix<- as.integer(output$total)
  fitted.ix<- as.integer(output$fitted)
  
  result <- cbind(knee.ix, inflection.ix, rank.ix, total.ix, fitted.ix)
  colnames(result) <- c("dropletUtils_barcodeRank_knee",
                        "dropletUtils_barcodeRank_inflection",
                        "dropletUtils_barcodeRank_rank",
                        "dropletUtils_barcodeRank_total",
                        "dropletUtils_barcodeRank_fitted")
  result.list <- list(result,
                      S4Vectors::metadata(output)$knee,
                      S4Vectors::metadata(output)$inflection)
  names(result.list) <- c("matrix","knee","inflection")
  return(result.list)
}


#' @title Identify empty droplets using \link[DropletUtils]{barcodeRanks}.
#' @description Run \link[DropletUtils]{barcodeRanks} on a count matrix
#'  provided in a \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain a raw counts matrix before empty droplets have been removed.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @param sample Character vector. Indicates which sample each cell belongs to
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample separately.
#'  If NULL, then all cells will be processed together. Default \code{NULL}.
#' @param lower See \link[DropletUtils]{emptyDrops} for more information. Default \code{100}.
#' @param fitBounds See \link[DropletUtils]{emptyDrops} for more information. Default \code{NULL}.
#' @param df See \link[DropletUtils]{emptyDrops} for more information. Default \code{20}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{barcodeRanks} output table appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{dropletUtils_BarcodeRank_Knee} and \emph{dropletUtils_BarcodeRank_Knee}
#'  Please refer to the documentation of \link[DropletUtils]{barcodeRanks} for
#'  details.
#' @examples
#' # The following unfiltered PBMC_1k_v3 data were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' # /pbmc_1k_v3
#' # Only the top 10 cells with most counts and the last 10 cells with non-zero
#' # counts are included in this example.
#' # This example only serves as an proof of concept and a tutoriol on how to
#' # run the function. The results should not be
#' # used for drawing scientific conclusions.
#' data(scExample, package = "singleCellTK")
#' sce <- runBarcodeRankDrops(inSCE = sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
runBarcodeRankDrops <- function(inSCE,
                                sample = NULL,
                                useAssay = "counts",
                                lower = 100,
                                fitBounds = NULL,
                                df = 20
) {
  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  } else {
    sample = rep(1, ncol(inSCE))
  }
  
  message(paste0(date(), " ... Running 'barcodeRanks'"))
  
  ##  Getting current arguments values
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  
  rank <- list()
  
  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
                                 dropletUtils_BarcodeRank_Knee = integer(ncol(inSCE)),
                                 dropletUtils_BarcodeRank_Inflection = integer(ncol(inSCE)))
  
  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  metaOutList <- list()
  for (i in seq_len(length(samples))) {
    sceSampleInd <- sample == samples[i]
    sceSample <- inSCE[, sceSampleInd]
    
    ## Define meta matrix for each subinSCE
    metaOutput <- S4Vectors::DataFrame(row.names = colnames(sceSample),
                                       dropletUtils_barcodeRank_rank = integer(ncol(sceSample)),
                                       dropletUtils_barcodeRank_total = integer(ncol(sceSample)),
                                       dropletUtils_barcodeRank_fitted = integer(ncol(sceSample)),
                                       dropletUtils_barcodeRank_knee = integer(ncol(sceSample)),
                                       dropletUtils_barcodeRank_inflection = integer(ncol(sceSample)),
                                       sample = colData(sceSample)[["sample"]])
    
    mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
    result <- .runBarcodeRankDrops(barcode.matrix = mat, lower=lower,
                                   fit.bounds=fitBounds,
                                   df=df)
    
    result.matrix <- result$matrix
    output[sceSampleInd, ] <- result.matrix[, c("dropletUtils_barcodeRank_knee", "dropletUtils_barcodeRank_inflection")]
    
    metaCols <- c("dropletUtils_barcodeRank_rank", "dropletUtils_barcodeRank_total", 
                  "dropletUtils_barcodeRank_fitted")
    metaOutput[sceSampleInd, metaCols] <- result.matrix[, metaCols]
    metaOutput[sceSampleInd,"dropletUtils_barcodeRank_knee"] <- rep(result$knee, sum(sceSampleInd))
    metaOutput[sceSampleInd,"dropletUtils_barcodeRank_inflection"] <- rep(result$inflection, sum(sceSampleInd))
    
    # Remove duplicated Rank
    metaOutput <- metaOutput[!duplicated(metaOutput$dropletUtils_barcodeRank_rank), ]
    
    metaOutList[[samples[i]]] <- metaOutput
  }
  
  colData(inSCE) = cbind(colData(inSCE), output)
  S4Vectors::metadata(inSCE)$runBarcodeRanksMetaOutput <- metaOutList
  
  inSCE@metadata$runBarcodeRankDrops <- argsList[-1]
  inSCE@metadata$runBarcodeRankDrops$packageVersion <- utils::packageDescription("DropletUtils")$Version
  
  return(inSCE)
}
