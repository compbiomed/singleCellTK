#' get Biomarker
#'
#' Use this function to get expression or binary data of gene list
#'
#' @param count_data A SCE object
#' @param gene gene list
#' @param binary binary or gradient
#'
#' @return A PCA plot
#' @export getBiomarker
#'

getBiomarker <- function(count_data, gene, binary="Binary"){
  # Get sample names
  sample <- rownames(pData(count_data))
  # Get counts for gene in sample
  c <- counts(count_data)[c(gene), ]
  # If color scale is "yes"/"no"
  if (binary == "Binary"){
    expression <- c > 0
  }
  # If color scale is a continuouse scale bar
  else if (binary == "Continuous"){
    expression <- exprs(count_data)[c(gene), ]
  }
  # Make data frame with sample, counts
  bio <- data.frame(sample, expression)
  colnames(bio) <- c("sample", gene)
  return(bio)
}
