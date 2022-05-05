#' Run EnrichR on SCE object
#' @details 
#' EnrichR works by querying the specified \code{features} to its online 
#' databases, thus it requires the Internet connection. 
#' 
#' Available \code{db} options could be shown by running 
#' \code{enrichR::listEnrichrDbs()$libraryName}
#' 
#' This function checks for the existence of features in the SCE object. When 
#' \code{features} do not have a match in \code{rownames(inSCE)}, users may 
#' try to specify \code{by} to pass the check. 
#' 
#' EnrichR expects gene symbols/names as the input (i.e. Ensembl ID might not 
#' work). When specified \code{features} are not qualified for this, users may 
#' try to specify \code{featureName} to change the identifier type to pass to 
#' EnrichR. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param features Character vector, selected genes for enrichment analysis. 
#' @param analysisName A string that identifies each specific analysis.
#' @param db Character vector. Selected database name(s) from the enrichR 
#' database list. If \code{NULL} then EnrichR will be run on all the available 
#' databases on the enrichR database. See details. Default \code{NULL}
#' @param by Character. From where should we find the \code{features}? 
#' \code{"rownames"} for from \code{rownames(inSCE)}, otherwise, from a column
#' of feature metadata (\code{rowData(inSCE)[[by]]}). See details. Default 
#' \code{"rownames"}.
#' @param featureName Character. Indicates the actual feature identifiers to be
#' passed to EnrichR. Can be \code{"rownames"}, a column in feature metadata 
#' (\code{rowData(inSCE)[[featureName]]}), or a character vector with its length
#' equals to \code{nrow(inSCE)}. See details. Default \code{"rownames"}.
#' @return Updates \code{inSCE} metadata with a data.frame of enrichment terms 
#' overlapping in the respective databases along with p-values, z-scores etc.
#' @export
#' @seealso \code{\link{getEnrichRResult}}
#' @examples
#' data("mouseBrainSubsetSCE")
#' if (Biobase::testBioCConnection()) {
#'   mouseBrainSubsetSCE <- runEnrichR(mouseBrainSubsetSCE, features = "Cmtm5", 
#'                                     db = "GO_Cellular_Component_2017",
#'                                     analysisName = "analysis1")
#' }
#' 
runEnrichR <- function(inSCE, 
                       features,
                       analysisName,
                       db = NULL,
                       by = "rownames",
                       featureName = NULL) {
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("inSCE has to inherit from SingleCellExperiment object.")
  }
  if (is.null(analysisName)) {
    stop("Have to specify analysisName.")
  }
  if (by == "rownames") {
    if (!all(features %in% rownames(inSCE))) {
      stop("Not all features found in `rownames(inSCE)`.")
    }
    allFeatures <- rownames(inSCE)
  } else {
    if (!by %in% names(SummarizedExperiment::rowData(inSCE))) {
      stop("`by` not found in rowData(inSCE).")
    }
    if (!all(features %in% SummarizedExperiment::rowData(inSCE)[[by]])) {
      stop("Not all features found in `rowData(inSCE)$",by,"`.")
    }
    allFeatures <- SummarizedExperiment::rowData(inSCE)[[by]]
  }
  if (!is.null(featureName)) {
    featureIdx <- allFeatures %in% features
    if (length(featureName) == 1) {
      if (featureName == "rownames") {
        features <- rownames(inSCE)[featureIdx]
      } else if (featureName %in% names(SummarizedExperiment::rowData(inSCE))) {
        features <- SummarizedExperiment::rowData(inSCE[featureIdx,])[[featureName]]
      } else {
        stop("featureName not found in `rowData(inSCE)`.")
      }
    } else if (length(featureName) == nrow(inSCE)) {
      features <- featureName[featureIdx]
    } else {
      stop("Invalid featureName specification.")
    }
  }
  
  internetConnection <- suppressWarnings(Biobase::testBioCConnection())
  #check for internet connection
  if (!internetConnection){
    stop("Please connect to the Internet and continue..")
  }
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  temp_db <- enrichR::listEnrichrDbs()
  enrdb <- temp_db$libraryName
  
  #test for db existing
  if (is.null(db)){
    db <- enrdb
  } else if (!all(db %in% enrdb)){
    db.notFound <- db[!db %in% enrdb]
    stop("database ", paste(db.notFound, collapse = ", "), " do not exist.")
  }
  
  enriched <- enrichR::enrichr(features, db)
  enriched <- data.frame(data.table::rbindlist(enriched, use.names = TRUE,
                                               fill = TRUE,
                                               idcol = "Database_selected"))
  
  enriched$link <- vapply(enriched$Database_selected, function(x){
    temp_db$link[temp_db$libraryName %in% x]
  }, FUN.VALUE = character(1))
  
  #sort the results based on p-values
  enriched <- enriched[order(enriched$P.value, decreasing = FALSE), ]
  
  #round the numeric values to their 7th digit
  #nums <- vapply(enriched, is.numeric, FUN.VALUE = logical(1))
  #enriched[, nums] <- round(enriched[, nums], digits = 7)
  
  getEnrichRResult(inSCE, analysisName) <- list(result = enriched,
                                                param = list(
                                                  features = features,
                                                  by = by,
                                                  db = db
                                                )) 
  return(inSCE)
}

#' @title Get or Set EnrichR Result
#' @rdname getEnrichRResult
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param analysisName A string that identifies each specific analysis
#' @param value The EnrichR result table
#' @return For getter method, a data.frame of the EnrichR result;
#' For setter method, \code{inSCE} with EnrichR results updated.
#' @export
#' @seealso \code{\link{runEnrichR}}
#' @examples 
#' data("mouseBrainSubsetSCE")
#' if (Biobase::testBioCConnection()) {
#'   mouseBrainSubsetSCE <- runEnrichR(mouseBrainSubsetSCE, features = "Cmtm5", 
#'                                     db = "GO_Cellular_Component_2017",
#'                                     analysisName = "analysis1")
#'   result <- getEnrichRResult(mouseBrainSubsetSCE, "analysis1")
#' }
setGeneric("getEnrichRResult<-", function(inSCE, analysisName, value) 
  standardGeneric("getEnrichRResult<-"))

#' @rdname getEnrichRResult
#' @export
setGeneric("getEnrichRResult", function(inSCE, analysisName) 
  standardGeneric("getEnrichRResult"))

#' @rdname getEnrichRResult
#' @export
setMethod("getEnrichRResult", 
          "SingleCellExperiment", 
          function(inSCE, analysisName){
            if (!"runEnrichR" %in% names(S4Vectors::metadata(inSCE)$sctk)) {
              stop("EnrichR analysis not performed yet. ",
                   "Please run `runEnrichR()`")
            }
            if (!analysisName %in% names(S4Vectors::metadata(inSCE)$sctk$runEnrichR)) {
              stop('"', analysisName, '" not found in EnrichR analysis names.')
            }
            return(S4Vectors::metadata(inSCE)$sctk$runEnrichR[[analysisName]])
          })

#' @rdname getEnrichRResult
#' @export
setReplaceMethod("getEnrichRResult", 
                 c("SingleCellExperiment"), 
                 function(inSCE, analysisName, value) {
                   if (is.null(analysisName)) {
                     stop("Have to specify analysisName.")
                   }
                   S4Vectors::metadata(inSCE)$sctk$runEnrichR[[analysisName]] <- value
                   return(inSCE)
                 })
