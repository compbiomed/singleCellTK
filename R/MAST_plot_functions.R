#' thresholdGenes
#'
#' @param SCEdata SingleCelltkExperiment object
#' @export
thresholdGenes <- function(SCEdata){
  # data preparation
  expres <- log2(assay(SCEdata, "counts") + 1)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)

  SCE_new <- SCE_new[which(freq(SCE_new) > 0), ]
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  return(thres)
}

#' MAST Violin
#'
#' Run MAST analysis on a SingleCelltkExperiment object.
#'
#' @param SCEdata SingleCelltkExperiment object
#' @param fcHurdleSig The filtered result from hurdle model
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#' @param variable Select the condition of interest
#'
#' @export
MASTviolin <- function(SCEdata, fcHurdleSig, samplesize = 49, threshP=FALSE,
                       variable){
  expres <- log2(assay(SCEdata, "counts") + 1)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  assays(SCE_new) <- list(thresh = thres$counts_threshold, tpm = assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:49]
  flat_dat <- as(SCE_new[entrez_to_plot, ], "data.table")
  if (threshP){
    yvalue <- "thresh"
  }
  else{
    yvalue <- "tpm"
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
#' Run MAST analysis on a SingleCelltkExperiment object.
#'
#' @param SCEdata SingleCelltkExperiment object
#' @param fcHurdleSig The filtered result from hurdle model
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#' @param variable Select the condition of interest
#'
#' @export
MASTregression <- function(SCEdata, fcHurdleSig, samplesize = 49, threshP=FALSE,
                           variable){
  count <- assay(SCEdata, "counts")
  expres <- log2(count + 1)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)
  SCE_new_lm <- MAST::FromMatrix(count, colData(SCEdata), fdata)
  cdr2 <- colSums(assay(SCE_new_lm) > 0)
  colData(SCE_new)$cngeneson <- scale(cdr2)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]

  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20, min_per_bin = 30)
  assays(SCE_new) <- list(thresh = thres$counts_threshold, tpm = assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:49]
  flat_dat <- as(SCE_new[entrez_to_plot, ], "data.table")

  if (threshP){
    yvalue <- "thresh"
  }
  else{
    yvalue <- "tpm"
  }

  ggbase <- ggplot2::ggplot(flat_dat, ggplot2::aes_string(x = variable,
                                                          y = yvalue,
                                                          color = variable)) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~primerid, scale = "free_y", ncol = 7)

  flat_dat[, lmPred := lm(thresh~cngeneson + condition)$fitted, key = primerid]

  regressionplot <- ggbase +
    ggplot2::aes(x = cngeneson) +
    ggplot2::geom_line(ggplot2::aes(y = lmPred), lty = 1) +
    ggplot2::xlab("Standardized Cellular Detection Rate")
  return(regressionplot)
}
