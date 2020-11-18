#' saveBiomarkerRes
#' Save biomarker gene information with a custom name when provided with diffex results.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param diffex results table saved from the differential expression analysis.
#' Required.
#' @param biomarkerName name of the biomarker result to be saved under in rowData().
#' Required.
#' @param ntop number of top N genes. Default is 25. Required
#' @param logFC logfold-change cutoff applied to save biomarker results. Optional
#' @param pVal adjusted p-value cutoff. Optional
#' @param method name of the diffex method used to generate the results. Options
#' are DESeq2, Limma and ANOVA. Required
#'
#' @return a new SCE object with the diffex result saved in the rowData using
#' the "biomarkerName"
#' @export
#'
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- scDiffEx(mouseBrainSubsetSCE,
#'                 useAssay = "logcounts",
#'                 condition = "level1class",
#'                 ntop = length(rownames(mouseBrainSubsetSCE)),
#'                 usesig = FALSE,
#'                 diffexmethod = "limma")
#'sceObj <- saveDiffExResults(mouseBrainSubsetSCE, res, "Limma_Level1class", "limma")
#'bioMarkRes <- saveBiomarkerRes(sceObj, res, "bioMarker", "limma", logFC = 4, pVal = 0.05)

saveBiomarkerRes <- function(inSCE, diffex, biomarkerName, method, ntop = 25, logFC = NULL, pVal = NULL) {
  if (!(class(inSCE) == "SingleCellExperiment" | class(inSCE) == "SCtkExperiment")) {
    stop("Please use a singleCellTK or a SCtkExperiment object")
  }
  if (!all(method %in% c("DESeq2", "limma", "ANOVA"))) {
    stop("Please provide a valid DiffEx method")
  }
  if (!is.null(biomarkerName)) {
    #check for special characaters
    if (grepl(" |-|_", biomarkerName)){
      biomarkerName <- gsub(" ", "_", biomarkerName)
    } else {
      if (grepl("[[:punct:]]", biomarkerName)) {
        stop("Please provide a name without special characters.")
      }
    }
  } else {
    stop("Please provide a name for the results.")
  }
  if (!is.data.frame(diffex)) {
    diffex <- data.frame(diffex)
  }
  bmRes <- NULL
  pvalIndex <- which(grepl("*padj*", colnames(diffex)))
  logFCIndex <- which(grepl("*log*", colnames(diffex)))
  if (length(logFCIndex) == 0) {
    logFCIndex <- 1
  }

  if (method == 'ANOVA' & !is.null(logFC)) {
    stop("logFC is not applicable for ANOVA")
  } else {
    if (!is.null(pVal) & !is.null(logFC)) {
      bmRes <- diffex[(diffex[, logFCIndex] <= logFC) & (diffex[, pvalIndex] <= pVal), ]
    } else if (!is.null(pVal)) {
      bmRes <- diffex[(diffex[, pvalIndex] <= pVal), ]
    } else if (!is.null(logFC)) {
      bmRes <- diffex[(diffex[, logFCIndex] <= logFC), ]
    }
  }

  if (is.null(bmRes)) {
    bmRes <- diffex
  }

  if (nrow(bmRes) < ntop) {
    bmRes <- bmRes[seq_len(nrow(bmRes)), ]
  } else {
    bmRes <- bmRes[seq_len(ntop), ]
  }
  rowData(inSCE)[, biomarkerName] <- ifelse(rownames(inSCE) %in% rownames(bmRes), 1, 0)
  return(inSCE)
}
