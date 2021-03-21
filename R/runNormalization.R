
#' Wrapper function to run any of the integrated normalization/transformation
#' methods in the singleCellTK.
#'
#' @param inSCE input
#' @param normalizationMethod normal
#' @param useAssay use
#' @param normAssayName npr
#' @param scale scale
#' @param seuratScaleFactor factor
#'
#' @return object
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