#' Run EnrichR on SCE object
#' @details 
#' The list of gene to be analyzed could be specified in two ways. 
#' Either by directly passing a character vector to \code{glist} and leave 
#' \code{geneSetCollectionName} and \code{geneSetName} as \code{NULL}, or by 
#' specifying the latter two arguments and leave \code{glist} as \code{NULL}.
#' 
#' Available \code{db} options could be shown by running 
#' \code{enrichR::listEnrichrDbs()$libraryName}
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param glist Character vector, selected genes for enrichment analysis using. 
#' Default \code{NULL}.
#' @param geneSetCollectionName Character. The name of an imported geneset 
#' collection. Default \code{NULL}.
#' @param geneSetName Character. The name of a geneset that should be found in 
#' the geneset collection specified by \code{geneSetCollectionName}. Default 
#' \code{NULL}. 
#' @param db Character vector. Selected database name(s) from the enrichR 
#' database list. If \code{NULL} then enrichR will be run on all the available 
#' databases on the enrichR database. Default \code{NULL}
#' @return Updates \code{inSCE} metadata with a data.frame of enrichment terms 
#' overlapping in the respective databases along with p-values, z-scores etc.
#' @export
#' @seealso \code{\link{getEnrichRResult}}
#' @examples
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- runEnrichR(mouseBrainSubsetSCE, glist = "Cmtm5", 
#'                                   db = "GO_Cellular_Component_2017")
runEnrichR <- function(inSCE, 
                       glist=NULL,
                       geneSetCollectionName=NULL, 
                       geneSetName=NULL, 
                       db = NULL) {
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("inSCE has to inherit from SingleCellExperiment object.")
  }
  if (is.null(glist)) {
    if (is.null(geneSetCollectionName) | is.null(geneSetName)) {
      stop("Either `glist` or {`geneSetCollectionName` and `geneSetName`} ",
           "should be specified.")
    }
    gsc <- .getGeneSetCollection(inSCE, geneSetCollectionName)
    if (!geneSetName %in% names(gsc)) {
      stop("Specified `geneSetName` not found in given `geneSetCollectionName`")
    }
    gs <- gsc[[geneSetName]]
    glist <- gs@geneIds
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
  
  enriched <- enrichR::enrichr(glist, db)
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
  
  getEnrichRResult(inSCE) <- enriched
  return(inSCE)
}

#' @title Get or Set EnrichR Result
#' @rdname getEnrichRResult
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param value The EnrichR result table
#' @return For getter method, a data.frame of the EnrichR result;
#' For setter method, \code{inSCE} with EnrichR results updated.
#' @export
#' @seealso \code{\link{runEnrichR}}
#' @examples 
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- runEnrichR(mouseBrainSubsetSCE, glist = "Cmtm5", 
#'                                   db = "GO_Cellular_Component_2017")
#' result <- getEnrichRResult(mouseBrainSubsetSCE)
setGeneric("getEnrichRResult<-", function(inSCE, value) 
  standardGeneric("getEnrichRResult<-") )

#' @rdname getEnrichRResult
#' @export
setGeneric("getEnrichRResult", function(inSCE) 
  standardGeneric("getEnrichRResult") )

#' @rdname getEnrichRResult
#' @export
setMethod("getEnrichRResult", 
          "SingleCellExperiment", 
          function(inSCE){
            return(S4Vectors::metadata(sce)$sctk$runEnrichR)
          })

#' @rdname getEnrichRResult
#' @export
setReplaceMethod("getEnrichRResult", 
                 c("SingleCellExperiment"), 
                 function(inSCE, value) {
                   S4Vectors::metadata(sce)$sctk$runEnrichR <- value
                   return(inSCE)
                 })