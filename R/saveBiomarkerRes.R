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
#'bioMarkRes <- saveBiomarkerRes(sceObj, res, "bioMarker", "limma", logFC = 0, pVal = 0.05)

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
  if (method == "limma"){
    if (!is.null(pVal) & !is.null(logFC)) {
      bmRes <- diffex[which((diffex[, 1] <= logFC) & (diffex[, 5] <= pVal)), ]
    } else if (!is.null(pVal)) {
      bmRes <- diffex[(diffex[, 5] <= pVal), ][seq_len(ntop), ]
    } else if (!is.null(logFC)) {
      bmRes <- diffex[(diffex[, 1] <= logFC), ][seq_len(ntop), ]
    } else {
      bmRes <- diffex[seq_len(ntop), ]
    }
  } else if (method == "DESeq2") {
    if (!is.null(pVal) & !is.null(logFC)) {
      bmRes <- diffex[which((diffex[, 2] <= logFC) & (diffex[, 6] <= pVal)), ][seq_len(ntop), ]
    } else if (!is.null(pVal)) {
      bmRes <- diffex[(diffex[, 6] <= pVal), ][seq_len(ntop), ]
    } else if (!is.null(logFC)) {
      bmRes <- diffex[(diffex[, 2] <= logFC), ][seq_len(ntop), ]
    } else {
      bmRes <- diffex[seq_len(ntop), ]
    }
  } else {
    if (!is.null(pVal) & !is.null(logFC)) {
      #bmRes <- diffex[which((diffex[, 2] <= logFC) & (diffex[, 6] <= pVal)), ][seq_len(ntop), ]
      stop("logFC is not applicable for ANOVA")
    } else if (!is.null(pVal)) {
      bmRes <- diffex[(diffex[, 2] <= pVal), ][seq_len(ntop), ]
    } else if (!is.null(logFC)) {
      #bmRes <- diffex[(diffex[, 2] <= logFC), ][seq_len(ntop), ]
      stop("logFC is not applicable for ANOVA")
    } else {
      bmRes <- diffex[seq_len(ntop), ]
    }
  }
  rowData(inSCE)[, biomarkerName] <- ifelse(rownames(inSCE) %in% rownames(bmRes), 1, 0)
  return(inSCE)
}
