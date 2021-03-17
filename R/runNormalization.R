
#' Title
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
                             normalizationMethod,
                             useAssay,
                             normAssayName,
                             scale = FALSE,
                             seuratScaleFactor = 10,
                             log = FALSE,
                             log1p = FALSE,
                             sqrt = FALSE,
                             pseudocountsNorm = NULL,
                             pseudocountsTransform = NULL,
                             trim = NULL
                             ){
  seuratMethods <- c("LogNormalize", "CLR", "RC", "SCTransform")
  scaterMethods <- c("logNormCounts", "CPM")
  tempSCE <- NULL
  
  #Perform 'Normalization'
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
  }
  else{
    stop("Specified normalization method '", normalizationMethod, "' not found.")
  }
  
  #Perform 'Scale'
  if(scale){
    scaledAssay <- computeZScore(counts = assay(tempSCE, normAssayName))
  }
  
  #Perform 'Transformation'
  if(log){
    selectedAssay <- log2(selectedAssay)
  }
  
  if(log1p){
    selectedAssay <- log1p(selectedAssay)
  }
  
  if(sqrt){
    selectedAssay <- sqrt(selectedAssay)
  }
  
  #Perform 'Pseudocounts' before Normalization
  if(!is.null(pseudocountsNorm)){
    selectedAssay <- selectedAssay + pseudocountsNorm
  }
  
  #Perform 'Pseudocounts' before Transformation
  if(!is.null(pseudocountsTransform)){
    selectedAssay <- selectedAssay + pseudocountsTransform
  }
  
  #Perform 'Trim'
  if(!is.null(trim)){
    selectedAssay <- trimCounts(selectedAssay, c(input$trimUpperValueAssay, input$trimLowerValueAssay))
  }
  
  return(inSCE)
}