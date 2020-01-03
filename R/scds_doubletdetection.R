#' @title annotates doublets using scds
#' @description Single Cell Doublet Scoring: In-silico doublet annotation for single cell RNA sequencing data
#' @param sce SingleCellExperiment object. Must contain a raw counts matrix.
#' @param rseed a random number seed for the scds method
#' @return SingleCellExperiment object with the scds output stored in coldata: cxds_score, bcds_score, hybrid_score
#' @examples
#' sce <- runScdsDoubletAnnotation(sce)
#' @export
#' @import scds
runScdsDoubletAnnotation <- function(sce,rseed = 12345) 
{
  ## install missing packages
  if("scds" %in% rownames(installed.packages()) == FALSE) {install.packages("scds")}
  library(scds)
  set.seed(rseed)
  sce = cxds(sce)
  sce = bcds(sce)
  sce = cxds_bcds_hybrid(sce)
 
  return(sce)
}

###### Test ##### 
test_scds_doublets <- function() {
  library('scds')
  data("sce_chcl")
  sce = sce_chcl
  sc <- runScdsDoubletAnnotation(sce)
  head(cbind(colData(sc)$cxds_score,colData(sc)$bcds_score,colData(sc)$hybrid_score))
}







