#' @describeIn MAST Identify adaptive thresholds
#'
#' @return thresholdGenes(): list of thresholded counts (on natural scale),
#' thresholds, bins, densities estimated on each bin, and the original data from
#' MAST::thresholdSCRNACountMatrix
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- thresholdGenes(mouseBrainSubsetSCE)
#'
thresholdGenes <- function(inSCE, useAssay="logcounts"){
  # data preparation
  expres <- SummarizedExperiment::assay(inSCE, useAssay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                              fdata)

  SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCENew),
                                           nbins = 20, min_per_bin = 30)
  return(thres)
}

#' @describeIn MAST Visualize MAST results using violin plots
#'
#' @param fcHurdleSig The filtered result from hurdle model
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#'
#' @return MASTviolin(): A ggplot object of MAST violin plots.
#' @export
#'
MASTviolin <- function(inSCE, useAssay="logcounts", fcHurdleSig,
                       samplesize = 49, threshP=FALSE, condition){
  expres <- SummarizedExperiment::assay(inSCE, useAssay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                             fdata)
  SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCENew),
                                           nbins = 20, min_per_bin = 30)
  names(SummarizedExperiment::assays(SCENew))[1] <- useAssay
  SummarizedExperiment::assay(SCENew, "thresh") <- thres$counts_threshold
  entrezToPlot <- fcHurdleSig$Gene[seq_len(min(nrow(fcHurdleSig), samplesize))]
  flatDat <- methods::as(SCENew[entrezToPlot, ], "data.table")
  if (threshP){
    yvalue <- "thresh"
  }
  else{
    yvalue <- useAssay
  }
  violinplot <- ggplot2::ggplot(flatDat,
                                ggplot2::aes_string(x = condition, y = yvalue,
                                                    color = condition)) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y", ncol = 7) +
    ggplot2::geom_violin() +
    ggplot2::ggtitle("Violin Plot")
  return(violinplot)
}

#' @describeIn MAST Visualize MAST results using linear model plots
#'
#' @return MASTregression(): A ggplot object of MAST linear regression plots.
#' @export
#'
MASTregression <- function(inSCE, useAssay="logcounts", fcHurdleSig,
                           samplesize = 49, threshP=FALSE, condition){
  expres <- SummarizedExperiment::assay(inSCE, useAssay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCENew <- MAST::FromMatrix(expres, SingleCellExperiment::colData(inSCE),
                             fdata)
  cdr2 <- colSums(SummarizedExperiment::assay(SCENew) > 0)
  SummarizedExperiment::colData(SCENew)$cngeneson <- scale(cdr2)
  SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]

  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCENew),
                                           nbins = 20, min_per_bin = 30)
  names(SummarizedExperiment::assays(SCENew))[1] <- useAssay
  SummarizedExperiment::assay(SCENew, "thresh") <- thres$counts_threshold
  entrezToPlot <- fcHurdleSig$Gene[seq_len(min(nrow(fcHurdleSig), samplesize))]

  flatDat <- methods::as(SCENew[entrezToPlot, ], "data.table")

  if (threshP){
    yvalue <- "thresh"
  } else{
    yvalue <- useAssay
  }

  resData <- NULL
  for (i in unique(flatDat$primerid)){
    resdf <- flatDat[flatDat$primerid == i, ]
    resdf$lmPred <- stats::lm(
      stats::as.formula(paste0(yvalue, "~cngeneson+", condition)),
      data = flatDat[flatDat$primerid == i, ])$fitted
    if (is.null(resData)){
      resData <- resdf
    } else {
      resData <- rbind(resData, resdf)
    }
  }

  ggbase <- ggplot2::ggplot(resData, ggplot2::aes_string(x = condition,
                                                          y = yvalue,
                                                          color = condition)) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y", ncol = 7)
  regressionplot <- ggbase +
    ggplot2::aes_string(x = "cngeneson") +
    ggplot2::geom_line(ggplot2::aes_string(y = "lmPred"), lty = 1) +
    ggplot2::xlab("Standardized Cellular Detection Rate")
  return(regressionplot)
}
