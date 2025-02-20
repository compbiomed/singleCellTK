.runBarcodeRankDrops <- function(barcode.matrix, lower = lower,
                                 fit.bounds = fit.bounds,
                                 df = df) {
  
  ## Convert to sparse matrix if not already in that format
  barcode.matrix <- .convertToMatrix(barcode.matrix)
  
  output <- DropletUtils::barcodeRanks(m = barcode.matrix, lower = lower,
                                       fit.bounds = fit.bounds,
                                       df = df)
  
  knee.ix <- as.integer(output$total >= 
                          S4Vectors::metadata(output)$knee)
  inflection.ix <- as.integer(output$total >= 
                                S4Vectors::metadata(output)$inflection)

  rank.ix <- as.integer(output$rank)
  total.ix <- as.integer(output$total)
  
  result <- cbind(knee.ix, inflection.ix, rank.ix, total.ix)
  colnames(result) <- c("dropletUtils_barcodeRank_knee",
                        "dropletUtils_barcodeRank_inflection",
                        "dropletUtils_barcodeRank_rank",
                        "dropletUtils_barcodeRank_total")
  
  result.list <- list(result,
                      S4Vectors::metadata(output)$knee,
                      S4Vectors::metadata(output)$inflection)
  names(result.list) <- c("matrix", "knee", "inflection")
  return(result.list)
}


#' @title Identify empty droplets using \link[DropletUtils]{barcodeRanks}.
#' @description Run \link[DropletUtils]{barcodeRanks} on a count matrix
#' provided in a \linkS4class{SingleCellExperiment} object. Distinguish between 
#' droplets containing cells and ambient RNA in a droplet-based single-cell RNA 
#' sequencing experiment.
#' @param inSCE A \linkS4class{SingleCellExperiment} object. Must contain a raw 
#' counts matrix before empty droplets have been removed.
#' @param sample Character vector or colData variable name. Indicates which 
#' sample each cell belongs to. Default \code{NULL}.
#' @param useAssay A string specifying which assay in the SCE to use. Default 
#' \code{"counts"}
#' @param lower See \link[DropletUtils]{barcodeRanks} for more information. 
#' Default \code{100}.
#' @param fitBounds See \link[DropletUtils]{barcodeRanks} for more information. 
#' Default \code{NULL}.
#' @param df See \link[DropletUtils]{barcodeRanks} for more information. Default 
#' \code{20}.
#' @return A \linkS4class{SingleCellExperiment} object with the
#' \link[DropletUtils]{barcodeRanks} output table appended to the
#' \link{colData} slot. The columns include
#' \code{dropletUtils_BarcodeRank_Knee} and 
#' \code{dropletUtils_barcodeRank_inflection}. Please refer to the documentation
#' of \link[DropletUtils]{barcodeRanks} for details.
#' @seealso \code{\link[DropletUtils]{barcodeRanks}}, 
#' \code{\link{runDropletQC}}, \code{\link{plotBarcodeRankDropsResults}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- runBarcodeRankDrops(inSCE = sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<- assay
runBarcodeRankDrops <- function(inSCE,
                                sample = NULL,
                                useAssay = "counts",
                                lower = 100,
                                fitBounds = NULL,
                                df = 20
) {
  
  p <- paste0(date(), " ... Running 'barcodeRanks'")
  message(p)
  
  ##  Getting current arguments values
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  argsList <- argsList[!names(argsList) %in% c("inSCE")]
  argsList$packageVersion <- utils::packageDescription("DropletUtils")$Version
  
  sample <- .manageCellVar(inSCE, var = sample)
  if (is.null(sample)) {
    sample <- rep(1, ncol(inSCE))
  }
  
  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(
    row.names = colnames(inSCE),
    dropletUtils_BarcodeRank_Knee = integer(ncol(inSCE)),
    dropletUtils_BarcodeRank_Inflection = integer(ncol(inSCE))
  )
  
  ## Loop through each sample and run barcodeRank
  samples <- unique(sample)
  for (s in samples) {
    sceSampleInd <- sample == s
    sceSample <- inSCE[, sceSampleInd]
    
    ## Define meta matrix for each subinSCE
    metaOutput <- S4Vectors::DataFrame(
      row.names = colnames(sceSample),
      dropletUtils_barcodeRank_rank = integer(ncol(sceSample)),
      dropletUtils_barcodeRank_total = integer(ncol(sceSample)),
      dropletUtils_barcodeRank_fitted = integer(ncol(sceSample)),
      dropletUtils_barcodeRank_knee = integer(ncol(sceSample)),
      dropletUtils_barcodeRank_inflection = integer(ncol(sceSample))
    )
    metaOutput$sample <- colData(sceSample)[["Sample"]]
    
    mat <- assay(sceSample, i = useAssay)
    result <- .runBarcodeRankDrops(barcode.matrix = mat, lower = lower,
                                   fit.bounds = fitBounds,
                                   df = df)
    
    result.matrix <- result$matrix
    output[sceSampleInd, ] <- 
      result.matrix[, c("dropletUtils_barcodeRank_knee",
                        "dropletUtils_barcodeRank_inflection")]
    
    metaCols <- c("dropletUtils_barcodeRank_rank", 
                  "dropletUtils_barcodeRank_total")
    metaOutput[, metaCols] <- result.matrix[, metaCols]
    metaOutput[,"dropletUtils_barcodeRank_knee"] <- rep(result$knee,
                                                        sum(sceSampleInd))
    metaOutput[,"dropletUtils_barcodeRank_inflection"] <- rep(result$inflection,
                                                              sum(sceSampleInd))
    
    # Remove duplicated Rank
    metaOutput <- 
      metaOutput[!duplicated(metaOutput$dropletUtils_barcodeRank_rank), ]
    if (!identical(samples, 1)) {
      S4Vectors::metadata(inSCE)$sctk$runBarcodeRankDrops[[s]] <- 
        list(metaOutput = metaOutput, argsList = argsList)
    }
  }
  if (identical(samples, 1)) {
    S4Vectors::metadata(inSCE)$sctk$runBarcodeRankDrops$all_cells <- 
      list(metaOutput = metaOutput, argsList = argsList)
  }
  
  colData(inSCE) <- cbind(colData(inSCE), output)
  return(inSCE)
}
