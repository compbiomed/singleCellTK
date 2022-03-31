#' @title Run VAM to score gene sets in single cell data
#' @description Wrapper for the Variance-adjusted Mahalanobis (VAM), which is a 
#' fast and accurate method for cell-specific gene set scoring of single cell 
#' data. This algorithm computes distance statistics and one-sided p-values for 
#' all cells in the specified single cell gene expression matrix. Gene sets 
#' should already be imported and stored in the meta data using functions such 
#' as \link{importGeneSetsFromList} or \link{importGeneSetsFromMSigDB}
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param geneSetCollectionName Character. The name of the gene set collection 
#' to use.
#' @param useAssay Character. The name of the assay to use. This assay should 
#' contain log normalized counts. 
#' @param resultNamePrefix  Character. Prefix to the name the VAM results which 
#' will be stored in the reducedDim slot of \code{inSCE}. The names of the 
#' output matrices will be \code{resultNamePrefix_Distance} and 
#' \code{resultNamePrefix_CDF}. If this parameter is set to \code{NULL}, then 
#' \code{"VAM_geneSetCollectionName_"} will be used. Default \code{NULL}.
#' @param center Boolean. If \code{TRUE}, values will be mean centered when 
#' computing the Mahalanobis statistic. Default \code{TRUE}.
#' @param gamma Boolean. If \code{TRUE}, a gamma distribution will be fit to 
#' the non-zero squared Mahalanobis distances computed from a row-permuted 
#' version of the gene expression matrix. The estimated gamma distribution will 
#' be used to compute a one-sided p-value for each cell. If \code{FALSE}, the 
#' p-value will be computed using the standard chi-square approximation for the 
#' squared Mahalanobis distance (or non-central if \code{center = FALSE}). 
#' Default \code{FALSE}.
#' @importFrom methods slot
#' @return A \linkS4class{SingleCellExperiment} object with VAM metrics stored 
#' in \code{reducedDim} as \code{VAM_NameOfTheGeneset_Distance} and 
#' \code{VAM_NameOfTheGeneset_CDF}.
#' @seealso \link{importGeneSetsFromList}, \link{importGeneSetsFromMSigDB}, 
#' \link{importGeneSetsFromGMT}, \link{importGeneSetsFromCollection} for 
#' importing gene sets. \link{sctkListGeneSetCollections}, 
#' \link{getPathwayResultNames} and \link{getGenesetNamesFromCollection} for 
#' available related information in \code{inSCE}.
#' @author Nida Pervaiz
#' @export
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' gs1 <- rownames(sce)[seq(10)]
#' gs2 <- rownames(sce)[seq(11,20)]
#' gs <- list("geneset1" = gs1, "geneset2" = gs2)
#' sce <- importGeneSetsFromList(inSCE = sce,geneSetList = gs,
#'                               by = "rownames")
#' sce <- runVAM(inSCE = sce, 
#'               geneSetCollectionName = "GeneSetCollection", 
#'               useAssay = "logcounts")
runVAM <- function(inSCE, geneSetCollectionName, useAssay, 
                   resultNamePrefix = NULL, center = TRUE, gamma = FALSE) { 
  ###################################################
  ###  create gene set collection 
  ###################################################
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("inSCE has to inherit from SingleCellExperiment object.")
  }
  gene.Set <- .getGeneSetCollection(inSCE, geneSetCollectionName)
  
  num.Genes <- length(gene.Set)
  gene.Set.Rows <- list()
  gene.Set.Collection <- list()
  
  for (i in seq(num.Genes)){
    
    gene.Set.Rows[i] <- slot(gene.Set[[i]], "setName")
    gene.Set.Ids <- slot(gene.Set[[i]], "geneIds")
    gene.Set.Collection[[i]] <- gene.Set.Ids
    
  }
  names(gene.Set.Collection) <- gene.Set.Rows
  
  gene.Set.Collection <- VAM::createGeneSetCollection(gene.ids = rownames(inSCE),
                                                      gene.set.collection = gene.Set.Collection)
  
  ###################################################
  ### execute vam 
  ###################################################
  mat <- expData(inSCE, useAssay)
  resultsexp <- VAM::vamForCollection(gene.expr = t(mat),
                                      gene.set.collection = gene.Set.Collection,
                                      center = center, gamma = gamma)
  
  ###################################################
  #### store results in metadata
  ###################################################
  if(is.null(resultNamePrefix)) {
    resultNamePrefix <- paste0("VAM_", geneSetCollectionName, "_")
  }
  SingleCellExperiment::reducedDim(inSCE, paste0(resultNamePrefix, "Distance")) <- resultsexp$distance.sq    
  SingleCellExperiment::reducedDim(inSCE, paste0(resultNamePrefix, "CDF")) <- resultsexp$cdf.value
  if ("pathwayAnalysisResultNames" %in% names(S4Vectors::metadata(inSCE))) {
    S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]] <- 
      c(S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]],
        paste0(resultNamePrefix, "Distance"),
        paste0(resultNamePrefix, "CDF"))
  } else {
    S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]] <- 
      c(paste0(resultNamePrefix, "Distance"),
        paste0(resultNamePrefix, "CDF"))
  }
  
  return (inSCE)
  
}
