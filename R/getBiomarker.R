#' Given a list of genes and a SingleCellExperiment object, return
#' the binary or continuous expression of the genes.
#'
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicates which assay to use. The default is "counts".
#' @param featureLocation Indicates which column name of rowData to query gene.
#' @param featureDisplay Indicates which column name of rowData to use
#' to display feature for visualization.
#'
#' @return getBiomarker(): A data.frame of expression values
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' getBiomarker(mouseBrainSubsetSCE, gene="C1qa")
#'
getBiomarker <- function(inSCE, gene, binary="Binary", useAssay="counts",
                         featureLocation = NULL, featureDisplay = NULL){
  # Get sample names
  sample <- colnames(inSCE)
  # Set rownames
  if(!is.null(featureLocation)){
    rownames(inSCE) = rowData(inSCE)[,featureLocation]
  }

  gene.ix <- c()
  for(g in gene){
    gene.ix <- c(gene.ix, which(rownames(inSCE) == g))
  }

  # Get counts for gene in sample
  c <- SummarizedExperiment::assay(inSCE, useAssay)[gene.ix, ,drop = FALSE]

  # If color scale is "yes"/"no"
  if (binary == "Binary"){
    expression <- c > 0
  } else if (binary == "Continuous"){
    expression <- c
  }
  # Make data frame with sample, counts
  bio <- cbind(sample, as.data.frame(t(as.matrix(expression))))

  if(!is.null(featureDisplay)){
    gene = rowData(inSCE)[gene.ix,featureDisplay]
  }
  colnames(bio) <- c("sample", gene)
  return(bio)
}
