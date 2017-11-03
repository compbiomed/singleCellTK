#' thresholdGenes
#'
#' @param SCEdata SCtkExperiment object
#' @param use_assay The assay to use for the MAST calculations. The default is
#' "logcounts"
#'
#' @export
thresholdGenes <- function(SCEdata, use_assay="logcounts"){
  # data preparation
  expres <- assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)

  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
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
#' @export
MASTviolin <- function(SCEdata, use_assay="logcounts", fcHurdleSig,
                       samplesize = 49, threshP=FALSE, variable){
  expres <- assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]
  thres <- MAST::thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  assays(SCE_new) <- list(thresh = thres$counts_threshold, tpm = assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:min(nrow(fcHurdleSig), samplesize)]
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
#' @export
MASTregression <- function(SCEdata, use_assay="logcounts", fcHurdleSig,
                           samplesize = 49, threshP=FALSE, variable){
  expres <- assay(SCEdata, use_assay)
  fdata <- data.frame(Gene = rownames(expres))
  rownames(fdata) <- fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, colData(SCEdata), fdata)
  cdr2 <- colSums(assay(SCE_new) > 0)
  colData(SCE_new)$cngeneson <- scale(cdr2)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new) > 0), ]

  thres <- MAST::thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                           min_per_bin = 30)
  assays(SCE_new) <- list(thresh = thres$counts_threshold, tpm = assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:min(nrow(fcHurdleSig), samplesize)]

  flat_dat <- as(SCE_new[entrez_to_plot, ], "data.table")

  if (threshP){
    yvalue <- "thresh"
  } else{
    yvalue <- "tpm"
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
