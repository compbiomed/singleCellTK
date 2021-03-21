
#' Wrapper function to run any of the integrated normalization/transformation
#' methods in the singleCellTK. The available methods include 'LogNormalize',
#' 'CLR', 'RC' and 'SCTransform' from Seurat, 'logNormCounts and 'CPM' from
#' Scater. Additionally, users can 'scale' using Z.Score, 'transform' using
#' log, log1p and sqrt, add 'pseudocounts' and trim the final matrices
#' between a range of values.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param normalizationMethod Specify a normalization method from `LogNormalize`,
#'  `CLR`, `RC` and `SCTransform` from Seurat or `logNormCounts` and `CPM` from
#'  scater packages. Default \code{NULL} is set which will not run any
#'  normalization method.
#' @param useAssay Specify the name of the assay that should be used.
#' @param normAssayName Specify the name of the new output assay.
#' @param scale Logical value indicating if the data should be scaled using
#'  Z.Score. Default \code{FALSE}.
#' @param seuratScaleFactor Specify the `scaleFactor` argument if a Seurat
#'  normalization method is selected. Default is \code{10}. This parameter
#'  will not be used if methods other than seurat are selected.
#' @param transformation Specify the transformation options to run on the
#'  selected assay. Options include `log2` (base 2 log transformation),
#'  `log1p` (natural log + 1 transformation) and `sqrt` (square root). Default
#'  value is \code{NULL}, which will not run any transformation. 
#' @param pseduoCountsNorm Specify a numeric pseudo value that should be added 
#'  to the assay before normalization is performed. Default is \code{NULL},
#'  which will not add any value.
#' @param pseduoCountsTransform Specify a numeric pseudo value that should be 
#'  added to the assay before transformation is run. Default is \code{NULL},
#'  which will not add any value.
#' @param trim Specify a vector of two numeric values that should be used
#'  as the upper and lower trim values to trim the assay between these two
#'  values. For example, \code{c(10,-10)} will trim the values between 10
#'  and -10. Default is \code{NULL}, which will not trim the data assay.
#'
#' @return Output SCE object with new normalized/transformed assay stored.
#' @export
#'
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- runNormalization(
#'  inSCE = sce_chcl, 
#'  "SCTransform", 
#'  "counts", 
#'  "sctCounts")
runNormalization <- function(inSCE,
                             normalizationMethod = NULL,
                             useAssay = "counts",
                             normAssayName = "customNormalizedAssay",
                             scale = FALSE,
                             seuratScaleFactor = 10,
                             transformation = NULL,
                             pseudocountsNorm = NULL,
                             pseudocountsTransform = NULL,
                             trim = NULL
                             ){
  seuratMethods <- c("LogNormalize", "CLR", "RC", "SCTransform")
  scaterMethods <- c("logNormCounts", "CPM")
  tempAssay <- NULL
  
  #Perform 'Pseudocounts' before Normalization
  if(!is.null(pseudocountsNorm)){
    tempAssay <- assay(inSCE, useAssay)
    tempAssay <- tempAssay + pseudocountsNorm
    assay(inSCE, useAssay) <- tempAssay
  }
  
  if(is.null(normalizationMethod)){
    #No normalizationMethod selected - Select the useAssay to perform other transformations
    tempAssay <- assay(inSCE, useAssay)
  }
  else{
    tempSCE <- NULL
    
    #Perform 'Normalization' - normalizationMethod is selected
    if(normalizationMethod %in% seuratMethods){
      if(normalizationMethod == "SCTransform"){
        tempSCE <- seuratSCTransform(
          inSCE = inSCE,
          normAssayName = normAssayName,
          useAssay = useAssay
        )
      }
      else{
        tempSCE <- seuratNormalizeData(
          inSCE = inSCE,
          normalizationMethod = normalizationMethod,
          useAssay = useAssay,
          normAssayName = normAssayName,
          scaleFactor = seuratScaleFactor
        )
      }
      tempAssay <- assay(tempSCE, normAssayName)
    }
    else if(normalizationMethod %in% scaterMethods){
      tempSCE <- do.call(
        paste0("scater_", normalizationMethod),
        list(
          inSCE = inSCE,
          assayName = normAssayName,
          useAssay = useAssay
        )
      )
      tempAssay <- assay(tempSCE, normAssayName)
    }
    else{
      stop("Specified normalization method '", normalizationMethod, "' not found.")
    }
  }
  
  #Perform 'Scale'
  if(scale){
    tempAssay <- computeZScore(counts = tempAssay)
  }
  
  #Perform 'Pseudocounts' before Transformation
  if(!is.null(pseudocountsTransform)){
    tempAssay <- tempAssay + pseudocountsTransform
  }
  
  #Perform 'Transformation'
  if("log2" %in% transformation){
    tempAssay <- log2(tempAssay)
  }
  
  if("log1p" %in% transformation){
    tempAssay <- log1p(tempAssay)
  }
  
  if("sqrt" %in% transformation){
    tempAssay <- sqrt(tempAssay)
  }
  
  #Perform 'Trim'
  if(!is.null(trim)){
    tempAssay <- trimCounts(tempAssay, c(trim[1], trim[2]))
  }
  
  assay(inSCE, normAssayName) <- tempAssay
  
  return(inSCE)
}