#' getBiomarker
#'
#' Given a list of genes and a SCtkExperiment object, return the binary or
#' continuous expression of the genes.
#'
#' @param count_data A SCtkExperiment object
#' @param gene gene list
#' @param binary "Binary" for binary expression or "Continuous" for a gradient.
#' Default: "Binary"
#'
#' @return A data.frame of expression values
#' @export getBiomarker
#'
getBiomarker <- function(count_data, gene, binary="Binary"){
  # Get sample names
  sample <- colnames(count_data)
  # Get counts for gene in sample
  c <- assay(count_data, "counts")[c(gene), ]
  # If color scale is "yes"/"no"
  if (binary == "Binary"){
    expression <- c > 0
  }
  # If color scale is a continuouse scale bar
  #TODO: change this to tpm assay or other normalized assay
  else if (binary == "Continuous"){
    expression <- log2(c + 1)
  }
  # Make data frame with sample, counts
  bio <- data.frame(sample, expression)
  colnames(bio) <- c("sample", gene)
  return(bio)
}
