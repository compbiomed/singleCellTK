#' getBiomarker
#'
#' Given a list of genes and a SCtkExperiment object, return the binary or
#' continuous expression of the genes.
#'
#' @param countData A SCtkExperiment object
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param useAssay Indicate which assay to use for PCA. Default is "logcounts"
#'
#' @return A data.frame of expression values
#' @export
#' @examples
#' getBiomarker(mouseBrainSubsetSCE, gene="C1qa")
#'
getBiomarker <- function(countData, gene, binary="Binary", useAssay="counts"){
  # Get sample names
  sample <- colnames(countData)
  # Get counts for gene in sample
  c <- SummarizedExperiment::assay(countData, useAssay)[c(gene), ]
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
