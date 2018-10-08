#' saveDiffExResults
#' Save Differential Expression Results with a custom name.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param diffex results table saved from the differential expression analysis.
#' Required
#' @param name name of the result to be saved under in rowData(). Required
#' @param method name of the diffex method used to generate the results. Options
#' are DESeq2, limma and ANOVA. Required
#'
#' @return a new SCE object with the diffex result saved in the rowData using the "name"
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffEx(mouseBrainSubsetSCE,
#'                 useAssay = "logcounts",
#'                 condition = "level1class",
#'                 ntop = length(rownames(mouseBrainSubsetSCE)),
#'                 usesig = FALSE,
#'                 diffexmethod = "limma")
#' sceObj <- saveDiffExResults(mouseBrainSubsetSCE, res, "Limma_Level1class", "limma")

saveDiffExResults <- function(inSCE, diffex, name, method) {
  if (!(class(inSCE) == "SingleCellExperiment" | class(inSCE) == "SCtkExperiment")) {
    stop("Please use a singleCellTK or a SCtkExperiment object")
  }
  if (!all(method %in% c("DESeq2", "limma", "ANOVA"))) {
    stop("Please provide a valid DiffEx method")
  }
  if (!is.null(name)) {
    #check for special characaters
    if (grepl(" |-|_", name)){
      name <- gsub(" ", "_", name)
    } else {
      if (grepl("[[:punct:]]", name)) {
        stop("Please provide a name without special characters.")
      }
    }
  } else {
    stop("Please provide a name for the results.")
  }
  #check if the results is in a data.frame() format
  if (!is.data.frame(diffex)) {
    diffex <- data.frame(diffex)
  }
  diffex <- diffex[rownames(inSCE), ]
  # if (any(rownames(diffex) != rownames(inSCE))) {
  #   stop("DiffEx Results and SCtkExperiment gene lists do not match.")
  # }
  colnames(diffex) <- paste((name), colnames(diffex), sep = "_")
  count_col <- ncol(rowData(inSCE))
  if (is.null(count_col)) {
    rowData(inSCE)[, seq_len(ncol(diffex))] <- diffex
  } else {
    rowData(inSCE)[, (count_col + seq_len(ncol(diffex)))] <- diffex
  }
  return(inSCE)
}
