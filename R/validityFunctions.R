.checkSCEValidity <- function(inSCE){
  if(is.null(rownames(inSCE)))
    stop("Rownames of the input SCE object cannot be NULL.")
  if(is.null(colnames(inSCE)))
    stop("Colnames of the input SCE object cannot be NULL.")
}

#' Check or retrieve cell metadata variable
#' @description Used as a helper for handling cell variables. Users can input
#' a custom variable or retrieve existing variable from \code{colData}. This
#' function makes sure that the output has valid length (\code{ncol(inSCE)}),
#' and convert vector to factor if needed.
#' @param inSCE Input \linkS4class{SingleCellExperiment} inherited object.
#' @param var A single character for \code{colData} variable or a vector of the
#' same length as \code{ncol(inSCE)}
#' @param as.factor Should the variable be output as \code{factor}?
#' @importFrom SummarizedExperiment colData
#' @return a vector of cell metadata variable, or a factor if
#' \code{as.factor = TRUE}.
#' @noRd
.manageCellVar <- function(inSCE, var = NULL, as.factor = FALSE) {
    if (!is.null(var)) {
      if (!is.vector(var) & !is.factor(var)) {
        stop("Invalid variable class")
      }
      if (length(var) != 1 & length(var) != ncol(inSCE)) {
        stop("Invalid variable length")
      }
      if (is.character(var)) {
        if (length(var) == 1) {
          if (!var %in% names(colData(inSCE))) {
            stop("Specified variable '", var, "' not found in colData(inSCE)")
          }
          var <- colData(inSCE)[[var]]
        }
      } else {
        if (length(var) == 1) {
          stop("Invalid variable class of length 1")
        }
      }
      if (isTRUE(as.factor) & !is.factor(var)) {
        var <- factor(var)
      }
    }
    return(var)
}

.manageFeatureVar <- function(inSCE, var = NULL, as.factor = FALSE) {
  if (!is.null(var)) {
    if (!is.vector(var) & !is.factor(var)) {
      stop("Invalid variable class")
    }
    if (length(var) != 1 & length(var) != nrow(inSCE)) {
      stop("Invalid variable length")
    }
    if (is.character(var)) {
      if (length(var) == 1) {
        if (!var %in% names(rowData(inSCE))) {
          stop("Specified variable '", var, "' not found in rowData(inSCE)")
        }
        var <- rowData(inSCE)[[var]]
      }
    } else {
      if (length(var) == 1) {
        stop("Invalid variable class of length 1")
      }
    }
    if (isTRUE(as.factor) & !is.factor(var)) {
      var <- factor(var)
    }
  }
  return(var)
}

#' Perform Checks and finally indicate which matrix to use.
#' @details Basic rule: (1) When \code{useAltExp} is not \code{NULL}, fetch
#' assay or reducedDim from altExp; (2) When \code{useReducedDim} is not
#' \code{NULL}, ignore \code{useAssay}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay Name of assay to use
#' @param useReducedDim Name of low-dimensional representation to use
#' @param useAltExp Name of alt-experiment to use
#' @param returnMatrix Whether to also return the matrix indicated by inputs.
#' @param cellAsCol When \code{returnMatrix = TRUE}, whether each column should
#' be a cell (i.e. transpose reducedDim).
#' @return A list, including \code{$names} for the final decision, and,
#' \code{$mat} for the matrix fetched when \code{returnMatrix = TRUE}.
#' @noRd
.selectSCEMatrix <- function(
    inSCE,
    useAssay = NULL,
    useReducedDim = NULL,
    useAltExp = NULL,
    returnMatrix = FALSE,
    cellAsCol = FALSE
) {
    if (!inherits(inSCE, "SingleCellExperiment")) {
        stop("`inSCE` should be a SingleCellExperiment Object.")
    }
    if (is.null(useAssay) & is.null(useReducedDim)) {
        stop("Either `useAssay` or `useReducedDim` has to be specified.")
    }
    if (is.null(useAltExp)) {
        if (!is.null(useReducedDim)) {
            if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)) {
                stop("Specified `useReducedDim` '", useReducedDim,
                     "' not found.")
            }
            useAssay <- NULL
        } else {
            if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)) {
                stop("Specified `useAssay` '", useAssay, "' not found.")
            }
        }
    } else {
        if (!useAltExp %in% SingleCellExperiment::altExpNames(inSCE)) {
            stop("Specified `useAltExp` '", useAltExp, "' not found.")
        }
        ae <- SingleCellExperiment::altExp(inSCE, useAltExp)
        if (!is.null(useReducedDim)) {
            if (!useReducedDim %in% SingleCellExperiment::reducedDimNames(ae)) {
                stop("Specified `useReducedDim` '", useReducedDim,
                     "' not found in the altExp.")
            }
            useAssay <- NULL
        } else {
            if(!useAssay %in% SummarizedExperiment::assayNames(ae)) {
                stop("Specified `useAssay` '", useAssay,
                     "' not found in the altExp.")
            }
        }
    }
    result <- list(
        names = list(
            useAssay = useAssay,
            useReducedDim = useReducedDim,
            useAltExp = useAltExp
        )
    )
    if (isTRUE(returnMatrix)) {
        if (!is.null(useAltExp)) {
            if (!is.null(useReducedDim)) {
                mat <- SingleCellExperiment::reducedDim(
                    SingleCellExperiment::altExp(inSCE, useAltExp),
                    useReducedDim)
                if (isTRUE(cellAsCol)) mat <- t(mat)
            } else {
                mat <- SummarizedExperiment::assay(
                    SingleCellExperiment::altExp(inSCE, useAltExp),
                    useAssay)
            }
        } else {
            if (!is.null(useReducedDim)) {
                mat <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)
                if (isTRUE(cellAsCol)) mat <- t(mat)
            } else {
                mat <- SummarizedExperiment::assay(inSCE, useAssay)
            }
        }
        result$mat <- mat
    }
    return(result)
}
