.myenv <- new.env(parent = emptyenv())

#' Run GSVA analysis on a \linkS4class{SingleCellExperiment} object
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param geneSetCollectionName Character. The name of the gene set collection to use.
#' parameter.
#' @param resultNamePrefix  Character. Prefix to the name the GSVA results which will be stored in the reducedDim slot of \code{inSCE}. The names of the output matrix will be \code{resultNamePrefix_Scores}. If this parameter is set to \code{NULL}, then "GSVA_geneSetCollectionName_" will be used. Default \code{NULL}.
#' @param ... Parameters to pass to gsva()
#'
#' @return   A \link[SingleCellExperiment]{SingleCellExperiment} object with pathway activity scores from GSVA stored in \code{reducedDim} as \code{GSVA_geneSetCollectionName_Scores}.

#' @export
#' 
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' gs1 <- rownames(sce)[seq(10)]
#' gs2 <- rownames(sce)[seq(11,20)]
#' gs <- list("geneset1" = gs1, "geneset2" = gs2)
#' 
#' sce <- importGeneSetsFromList(inSCE = sce,geneSetList = gs,
#'                                            by = "rownames")
#' sce <- runGSVA(inSCE = sce, geneSetCollectionName = "GeneSetCollection", useAssay = "logcounts")
#'
#'  

runGSVA <- function(inSCE, useAssay = "logcounts",
                    resultNamePrefix = NULL, geneSetCollectionName, ...){
    gsvaRes <- NULL
    gene.Set <- .getGeneSetCollection(inSCE, geneSetCollectionName)
  
  
    gsvaRes <- t(GSVA::gsva(as.matrix(expData(inSCE, useAssay)),
                                           gene.Set))
  
    
   if(is.null(resultNamePrefix)) {
      resultNamePrefix <- paste0("GSVA_", geneSetCollectionName, "_")
   }
    
   SingleCellExperiment::reducedDim(inSCE, paste0(resultNamePrefix, "Scores")) <- gsvaRes 
   
   if ("pathwayAnalysisResultNames" %in% names(S4Vectors::metadata(inSCE))) {
     S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]] <- 
       c(S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]],
         paste0(resultNamePrefix, "Scores"))
   } else {
     S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]] <- 
       c(paste0(resultNamePrefix, "Scores"))
   }
  
  
  return(inSCE)
}

