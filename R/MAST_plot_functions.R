#' thresholdGenes
#'
#' @param SCEdata SCtkExperiment object
#' @param use_assay The assay to use for the MAST calculations. The default is
#' "logcounts"
#'
#' @return list of thresholded counts (on natural scale), thresholds, bins,
#' densities estimated on each bin, and the original data from
#' MAST::thresholdSCRNACountMatrix
#' @export
#' @examples
#' data("GSE60361_subset_sce")
#' res <- thresholdGenes(GSE60361_subset_sce)
#'
thresholdGenes <- function(SCEdata, use_assay="logcounts"){
  # data preparation
  expres <- SummarizedExperiment::assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, SingleCellExperiment::colData(SCEdata),
                              fdata)

  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  return(thres)
}

#' MAST Violin
#'
#' Run MAST analysis on a SCtkExperiment object.
#'
#' @param SCEdata SCtkExperiment object
#' @param use_assay The assay to use for the MAST calculations. The default is
#' "logcounts"
#' @param fcHurdleSig The filtered result from hurdle model
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#' @param variable Select the condition of interest
#'
#' @return A ggplot object of MAST violin plots.
#' @export
#'
MASTviolin <- function(SCEdata, use_assay="logcounts", fcHurdleSig,
                       samplesize = 49, threshP=FALSE, variable){
  expres <- SummarizedExperiment::assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, SingleCellExperiment::colData(SCEdata), fdata)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  names(SummarizedExperiment::assays(SCE_new))[1] <- use_assay
  SummarizedExperiment::assay(SCE_new, "thresh") <- thres$counts_threshold
  entrez_to_plot <- fcHurdleSig$Gene[1:min(nrow(fcHurdleSig), samplesize)]
  flat_dat <- methods::as(SCE_new[entrez_to_plot, ], "data.table")
  if (threshP){
    yvalue <- "thresh"
  }
  else{
    yvalue <- use_assay
  }
  violinplot <- ggplot2::ggplot(flat_dat, ggplot2::aes_string(x = variable,
                                                              y = yvalue,
                                                              color = variable)) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y", ncol = 7) +
    ggplot2::geom_violin() +
    ggplot2::ggtitle("Violin Plot")
  return(violinplot)
}

#' MAST linear model
#'
#' Run MAST analysis on a SCtkExperiment object.
#'
#' @param SCEdata SCtkExperiment object
#' @param use_assay The assay to use for the MAST calculations. The default is
#' "logcounts"
#' @param fcHurdleSig The filtered result from hurdle model
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#' @param variable Select the condition of interest
#'
#' @return A ggplot object of MAST linear regression plots.
#' @export
#'
MASTregression <- function(SCEdata, use_assay="logcounts", fcHurdleSig,
                           samplesize = 49, threshP=FALSE, variable){
  expres <- SummarizedExperiment::assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, SingleCellExperiment::colData(SCEdata), fdata)
  cdr2 <- colSums(SummarizedExperiment::assay(SCE_new) > 0)
  SummarizedExperiment::colData(SCE_new)$cngeneson <- scale(cdr2)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]

  thres <- MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(SCE_new), nbins = 20,
                                           min_per_bin = 30)
  names(SummarizedExperiment::assays(SCE_new))[1] <- use_assay
  SummarizedExperiment::assay(SCE_new, "thresh") <- thres$counts_threshold
  entrez_to_plot <- fcHurdleSig$Gene[1:min(nrow(fcHurdleSig), samplesize)]

  flat_dat <- methods::as(SCE_new[entrez_to_plot, ], "data.table")

  if (threshP){
    yvalue <- "thresh"
  } else{
    yvalue <- use_assay
  }

  res_data <- NULL
  for (i in unique(flat_dat$primerid)){
    resdf <- flat_dat[flat_dat$primerid == i, ]
    resdf$lmPred <- lm(as.formula(paste0(yvalue, "~cngeneson+", variable)),
                       data = flat_dat[flat_dat$primerid == i, ])$fitted
    if (is.null(res_data)){
      res_data <- resdf
    } else {
      res_data <- rbind(res_data, resdf)
    }
  }

  ggbase <- ggplot2::ggplot(res_data, ggplot2::aes_string(x = variable,
                                                          y = yvalue,
                                                          color = variable)) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y", ncol = 7)
  regressionplot <- ggbase +
    ggplot2::aes(x = cngeneson) +
    ggplot2::geom_line(ggplot2::aes(y = lmPred), lty = 1) +
    ggplot2::xlab("Standardized Cellular Detection Rate")
  return(regressionplot)
}
