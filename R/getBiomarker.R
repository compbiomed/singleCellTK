#' Given a list of genes and a SingleCellExperiment object, return
#' the binary or continuous expression of the genes.
#'
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "counts".
#'
#' @return getBiomarker(): A data.frame of expression values
#' @export
#' @examples
#' getBiomarker(mouseBrainSubsetSCE, gene="C1qa")
#'
getBiomarker <- function(inSCE, gene, binary="Binary", useAssay="counts"){
  # Get sample names
  sample <- colnames(inSCE)
  # Get counts for gene in sample
  c <- SummarizedExperiment::assay(inSCE, useAssay)[c(gene), ]
  # If color scale is "yes"/"no"
  if (binary == "Binary"){
    expression <- c > 0
  } else if (binary == "Continuous"){
    expression <- c
  }
  # Make data frame with sample, counts
  bio <- data.frame(sample, expression)
  colnames(bio) <- c("sample", gene)
  return(bio)
}
