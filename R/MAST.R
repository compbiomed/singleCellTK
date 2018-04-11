#' MAST
#'
#' Run MAST analysis on a SCESet object.
#'
#' @param SCEdata SCESet object
#' @param useAssay The assay to use for the MAST calculations. The default is
#' "logcounts"
#' @param condition select variable (from the colData) that is used for the model.
#' @param interest.level If the condition of interest has more than two factors,
#' indicate which level should be used to compare to all other samples.
#' @param freqExpressed Filter genes that are expressed in at least this
#' fraction of cells. The default is expression in 0.1 of samples.
#' @param fcThreshold Minimum fold change for differentially expressed gene.
#' @param p.value p values for selecting the hurdle result, default is 0.05
#' @param useThresh Use adaptive thresholding to filter genes. The default is
#' FALSE.
#'
#' @return A data.frame of differentially expressed genes with p-values.
#' @export
#'
MAST <- function(SCEdata, condition = NULL, interest.level = NULL,
                 freqExpressed = 0.1, fcThreshold=log2(1.5), p.value = 0.05,
                 useThresh=FALSE, useAssay = "logcounts"){

  if (is.null(condition)){
    stop("specify the condition of interest")
  }

  if (length(unique(SingleCellExperiment::colData(SCEdata)[, condition])) == 1) {
    stop("only one level is in the condition")
  }

  if (is.null(interest.level) & length(
    unique(SingleCellExperiment::colData(SCEdata)[, condition])) > 2){
    stop("You must specify a level of interest when more than 2 levels are in",
         " the condition")
  }

  # Create MAST SingleCellAssay
  pdata <- SingleCellExperiment::colData(SCEdata)
  expres <- SummarizedExperiment::assay(SCEdata, useAssay)
  fdata <- SingleCellExperiment::rowData(SCEdata)
  SCENew <- MAST::FromMatrix(expres, pdata, fdata)

  #Caculate CDR for zlm model
  SummarizedExperiment::colData(SCENew)$cngeneson <-
    scale(colSums(SummarizedExperiment::assay(SCENew) > 0))

  if (useThresh){
    SCENew <- SCENew[which(MAST::freq(SCENew) > 0), ]
    thresh <- MAST::thresholdSCRNACountMatrix(
      SummarizedExperiment::assay(SCENew), nbins = 20, min_per_bin = 30)
    SummarizedExperiment::assays(SCENew) <-
      list(thresh = thresh$counts_threshold,
           tpm = SummarizedExperiment::assay(SCENew))
  }

  # filter based on frequency of expression across samples
  SCENewSample <- SCENew[which(MAST::freq(SCENew) > freqExpressed), ]

  # if the condition of interest is numeric, to change it to a factor
  if (is.numeric(SummarizedExperiment::colData(SCENewSample)[, condition])){
    SummarizedExperiment::colData(SCENewSample)[, condition] <-
      as.factor(SummarizedExperiment::colData(SCENewSample)[, condition])
  }

  # >2 levels in the condition
  if (!is.null(interest.level) &
      length(unique(SingleCellExperiment::colData(SCEdata)[, condition])) > 2){
    levels(SummarizedExperiment::colData(SCENewSample)[, condition]) <-
      c(levels(SummarizedExperiment::colData(SCENewSample)[, condition]),
        paste0("no_", interest.level))
    SummarizedExperiment::colData(SCENewSample)[, condition][
      SummarizedExperiment::colData(SCENewSample)[, condition] !=
        interest.level] <- paste0("no_", interest.level)
    SummarizedExperiment::colData(SCENewSample)[, condition] <-
      droplevels(as.factor(
        SummarizedExperiment::colData(SCENewSample)[, condition]))

    hurdle1 <- MAST::zlm(stats::as.formula(paste0("~", condition,
                                                  "+cngeneson")), SCENewSample)
    summaryh1 <- MAST::summary(hurdle1, doLRT = paste0(condition, "no_",
                                                       interest.level))

    summaryDT <- summaryh1[["datatable"]]

    fcHurdle <- merge(
      summaryDT[summaryDT$contrast == paste0(condition, "no_", interest.level) &
                  summaryDT$component == "H", c("primerid", "Pr(>Chisq)")],
      summaryDT[summaryDT$contrast == paste0(condition, "no_", interest.level) &
                  summaryDT$component == "logFC", c("primerid", "coef", "ci.hi",
                                                    "ci.lo")]
    )
  } else {
    SummarizedExperiment::colData(SCENewSample)[, condition] <-
      droplevels(as.factor(
        SummarizedExperiment::colData(SCENewSample)[, condition]))
    level.cond  <- levels(
      SummarizedExperiment::colData(SCENewSample)[, condition])

    hurdle1 <- MAST::zlm(stats::as.formula(paste0("~", condition,
                                                  "+cngeneson")), SCENewSample)
    summaryh1 <- MAST::summary(hurdle1, doLRT = paste0(condition,
                                                       level.cond[2]))

    summaryDT <- summaryh1[["datatable"]]

    fcHurdle <- merge(
      summaryDT[summaryDT$contrast == paste0(condition, level.cond[2]) &
                  summaryDT$component == "H", c("primerid", "Pr(>Chisq)")],
      summaryDT[summaryDT$contrast == paste0(condition, level.cond[2]) &
                  summaryDT$component == "logFC", c("primerid", "coef", "ci.hi",
                                                    "ci.lo")]
    )
  }

  # Use p-value correction method, here we use fdr
  fcHurdle$fdr <- stats::p.adjust(fcHurdle$"Pr(>Chisq)", "fdr")

  # Filter the data again by the adjusted pvalue and coef
  fcHurdleSig <-  fcHurdle[fcHurdle$fdr < p.value & abs(fcHurdle$coef) >
                             fcThreshold & !is.nan(fcHurdle$coef), ]
  colnames(fcHurdleSig)[1] <- "Gene"

  data.table::setorder(fcHurdleSig, "fdr")

  return(fcHurdleSig)
}
