
#' Wrapper function to run any of the integrated normalization/transformation
#' methods in the singleCellTK. The available methods include 'LogNormalize',
#' 'CLR', 'RC' and 'SCTransform' from Seurat, 'logNormCounts and 'CPM' from
#' Scater. Additionally, users can 'scale' using Z.Score, 'transform' using
#' log, log1p and sqrt, add 'pseudocounts' and trim the final matrices
#' between a range of values.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Specify the name of the assay that should be used.
#' @param outAssayName Specify the name of the new output assay.
#' @param normalizationMethod Specify a normalization method from `LogNormalize`,
#'  `CLR`, `RC` and `SCTransform` from Seurat or `logNormCounts` and `CPM` from
#'  scater packages. Default \code{NULL} is set which will not run any
#'  normalization method.
#' @param scale Logical value indicating if the data should be scaled using
#'  Z.Score. Default \code{FALSE}.
#' @param seuratScaleFactor Specify the `scaleFactor` argument if a Seurat
#'  normalization method is selected. Default is \code{10000}. This parameter
#'  will not be used if methods other than seurat are selected.
#' @param transformation Specify the transformation options to run on the
#'  selected assay. Options include `log2` (base 2 log transformation),
#'  `log1p` (natural log + 1 transformation) and `sqrt` (square root). Default
#'  value is \code{NULL}, which will not run any transformation.
#' @param pseudocountsBeforeNorm Specify a numeric pseudo value that should be added
#'  to the assay before normalization is performed. Default is \code{NULL},
#'  which will not add any value.
#' @param pseudocountsBeforeTransform Specify a numeric pseudo value that should be
#'  added to the assay before transformation is run. Default is \code{NULL},
#'  which will not add any value.
#' @param trim Specify a vector of two numeric values that should be used
#'  as the upper and lower trim values to trim the assay between these two
#'  values. For example, \code{c(10,-10)} will trim the values between 10
#'  and -10. Default is \code{NULL}, which will not trim the data assay.
#' @param verbose Logical value indicating if progress messages should be
#' displayed to the user. Default is \code{TRUE}.
#'
#' @return Output SCE object with new normalized/transformed assay stored.
#' @export
#'
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- runNormalization(
#'  inSCE = sce_chcl,
#'  normalizationMethod = "LogNormalize",
#'  useAssay = "counts",
#'  outAssayName = "logcounts")
runNormalization <- function(inSCE,
                             useAssay = "counts",
                             outAssayName = "customNormalizedAssay",
                             normalizationMethod = NULL,
                             scale = FALSE,
                             seuratScaleFactor = 10000,
                             transformation = NULL,
                             pseudocountsBeforeNorm = NULL,
                             pseudocountsBeforeTransform = NULL,
                             trim = NULL,
                             verbose = TRUE
                             ){
  seuratMethods <- c("LogNormalize", "CLR", "RC", "SCTransform")
  scaterMethods <- c("logNormCounts", "CPM")
  tempAssay <- NULL
  tag <- "normalized"

  #Perform 'Pseudocounts' before Normalization
  if(!is.null(pseudocountsBeforeNorm)){
    tempAssay <- assay(inSCE, useAssay)
    tempAssay <- tempAssay + pseudocountsBeforeNorm
    assay(inSCE, useAssay) <- tempAssay

    if(verbose)
      message("Added ", pseudocountsBeforeNorm, " to input matrix before normalizing data.")
  }

  #Normalization (if applicable)
  if(is.null(normalizationMethod)){
    #No normalizationMethod selected - Select the useAssay to perform other transformations
    tempAssay <- assay(inSCE, useAssay)
  }
  else{
    tempSCE <- NULL

    #Perform 'Normalization' - normalizationMethod is selected
    if(normalizationMethod %in% seuratMethods){
      if(normalizationMethod == "SCTransform"){
        tempSCE <- runSeuratSCTransform(
          inSCE = inSCE,
          normAssayName = outAssayName,
          useAssay = useAssay,
          verbose = verbose
        )
        tag <- "transformed"
      }
      else{
        tempSCE <- runSeuratNormalizeData(
          inSCE = inSCE,
          normalizationMethod = normalizationMethod,
          useAssay = useAssay,
          normAssayName = outAssayName,
          scaleFactor = seuratScaleFactor,
          verbose = verbose
        )
        if(normalizationMethod %in% c("LogNormalize", "CLR")){
          tag <- "transformed"
        }
      }
      tempAssay <- assay(tempSCE, outAssayName)
    }
    else if(normalizationMethod %in% scaterMethods){
      tempSCE <- do.call(
        paste0("scater", normalizationMethod),
        list(
          inSCE = inSCE,
          assayName = outAssayName,
          useAssay = useAssay
        )
      )
      tempAssay <- assay(tempSCE, outAssayName)
      if(normalizationMethod == "logNormCounts"){
        tag <- "transformed"
      }
    }
    else{
      stop("Specified normalization method '", normalizationMethod, "' not found.")
    }

    if(verbose)
      message("Normalization performed using ", normalizationMethod, " method.")
  }

  #Perform 'Pseudocounts' before Transformation
  if(!is.null(pseudocountsBeforeTransform)){
    tempAssay <- tempAssay + pseudocountsBeforeTransform

    if(verbose)
      message("Added ", pseudocountsBeforeTransform, " to input matrix before transforming data.")
  }

  #Perform 'Transformation'
  if("log2" %in% transformation){
    tempAssay <- log2(tempAssay)
    tag <- "transformed"

    if(verbose)
      message("Log2 transformation performed on the input data.")
  }

  if("log1p" %in% transformation){
    tempAssay <- log1p(tempAssay)
    tag <- "transformed"

    if(verbose)
      message("Log1p transformation (natural log + 1) performed on the input data.")
  }

  if("sqrt" %in% transformation){
    tempAssay <- sqrt(tempAssay)
    tag <- "transformed"

    if(verbose)
      message("Sqrt transformation (square root) performed on the input data.")
  }

  #Perform 'Scale'
  if(scale){
    tempAssay <- computeZScore(counts = tempAssay)
    tag <- "scaled"

    if(verbose)
      message("Z.Score scaling performed on the input data.")
  }

  #Perform 'Trim'
  if(!is.null(trim)){
    tempAssay <- trimCounts(tempAssay, c(trim[1], trim[2]))

    if(verbose)
      message("Data trimmed between ", trim[1], " and ", trim[2], ".")
  }

  expData(inSCE, outAssayName, tag = tag, altExp = FALSE) <- tempAssay

  return(inSCE)
}