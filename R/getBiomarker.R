#' getBiomarker
#'
#' Given a list of genes and a SCtkExperiment object, return the binary or
#' continuous expression of the genes.
#'
#' @param count_data A SCtkExperiment object
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#' @param use_assay Indicate which assay to use for PCA. Default is "logcounts"
#'
#' @return A data.frame of expression values
#' @export
#'
getBiomarker <- function(count_data, gene, binary="Binary", use_assay="counts"){
  # Get sample names
  sample <- colnames(count_data)
  # Get counts for gene in sample
  c <- assay(count_data, use_assay)[c(gene), ]
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
