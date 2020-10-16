#' @title Subset a SingleCellExperiment object by rows
#' @description Used to peform subsetting of a
#' \linkS4class{SingleCellExperiment} object using a variety of methods that
#' indicate the correct rows to keep. The various methods,
#' \code{index}, \code{bool}, and \code{rowData}, can be used in conjunction
#' with one another. If \code{returnAsAltExp} is set to \code{TRUE},
#' then the returned object will have the same number of rows as the input
#' \code{inSCE} as the subsetted object will be stored in the
#' \code{\link[SingleCellExperiment]{altExp}} slot.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param index Integer vector. Vector of indicies indicating which rows
#' to keep. If \code{NULL}, this will not be used for subsetting.
#' Default \code{NULL}.
#' @param bool Boolean vector. Vector of \code{TRUE} or \code{FALSE}
#' indicating which rows should be kept. Needs to be the same length as the
#' number of rows in \code{inSCE}. If \code{NULL}, this will not be used
#' for subsetting. Default \code{NULL}.
#' @param rowData Character. An expression that will identify a subset of rows
#' using variables found in the \code{rowData} of \code{inSCE}. For example,
#' if \code{x} is a numeric vector in \code{rowData}, then \code{"x < 5"} will
#' return all rows with x less than 5. Single quotes should be used for
#' character strings. For example, \code{"y == 'yes'"} will return all
#' rows where y is "yes". Multiple expressions can be evaluated by placing them
#' in a vector. For example \code{c("x < 5", "y =='yes'")} will apply both
#' operations for subsetting. If \code{NULL}, this will not be used for
#' subsetting. Default \code{NULL}.
#' @param returnAsAltExp Boolean. If \code{TRUE}, the subsetted
#' \linkS4class{SingleCellExperiment} object will be returned in the
#' \code{altExp} slot of \code{inSCE}. If \code{FALSE}, the subsetted
#' \linkS4class{SingleCellExperiment} object will be directly returned.
#' @param altExpName Character. Name of the alternative experiment object to
#' add if \code{returnAsAltExp = TRUE}. Default \code{subset}.
#' @param prependAltExpName Boolean. If \code{TRUE}, \code{altExpName} will
#' be added to the beginning of the assay names in the \code{altExp} object.
#' This is only utilized if \code{returnAsAltExp = TRUE}. Default \code{TRUE}.
#' @author Joshua D. Campbell
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object that has
#' been subsetted by rowData.
#' @examples
#' data(scExample)
#'
#' # Set a variable up in the rowData indicating mitochondrial genes
#' rowData(sce)$isMito <- ifelse(grepl("^MT-", rowData(sce)$feature_name),
#'                               "yes", "no")
#' sce <- subsetSCERows(sce, rowData = "isMito == 'yes'")
#' @export
#' @importFrom SummarizedExperiment assays assays<-
subsetSCERows <- function(inSCE, index = NULL, bool = NULL, rowData = NULL,
                          returnAsAltExp = TRUE, altExpName = "subset",
                          prependAltExpName = TRUE) {

  if(is.null(index) & is.null(bool) & is.null(rowData)) {
    stop("At least one of 'index', 'bool', or 'rowData' must be supplied.")
  }
  final.ix <- rep(FALSE, nrow(inSCE))

  # Parse index containing integers
  if(!is.null(index)) {
    if(min(index) < 1 | max(index) > nrow(inSCE)) {
      stop("'index' must contain integers between 1 and the number of rows ",
           "in 'inSCE': ", nrow(inSCE))
    }
    final.ix[index] <- TRUE
  }

  # Parse Boolean vector
  if(!is.null(bool)) {
    if(length(bool) != nrow(inSCE) | !is.logical(bool)) {
      stop("'bool' must be a logical vector the same length as the number of ",
           "rows in 'inSCE': ", nrow(inSCE))
    }
    final.ix[bool] <- TRUE
  }

  # Parse expressions for rowData variables
  if(!is.null(rowData)) {
    for(i in seq_along(rowData)) {
      temp <- eval(parse(text = as.character(rowData[i])),
                   envir = as.data.frame(rowData(inSCE)))
      if(length(temp) != nrow(inSCE) | !is.logical(temp)) {
        stop("The expression ", rowData[i], " did not produce a boolean ",
             "vector the same length as the number of rows in 'inSCE'. ",
             "Please ensure that the spelling of the variable you are ",
             "trying to use matches one of the column names in ",
             "'rowData(inSCE)' and that your expression is valid.")
      }
      final.ix[temp] = TRUE
    }
  }

  # Filter object and return in different ways depending on flag
  temp.SCE <- inSCE[final.ix,]
  if(isTRUE(returnAsAltExp)) {
    if(is.null(altExpName) || length(altExpName) > 1
                           || !is.character(altExpName)) {
      stop("'altExpName' needs to be a character of length 1 if the subset
           is going to be saved as an ", "'altExp' object.")
    }
    if(isTRUE(prependAltExpName)) {
      names(assays(temp.SCE)) <- paste0(altExpName, names(assays(temp.SCE)))
    }
    SingleCellExperiment::altExp(inSCE, altExpName) <- temp.SCE
  } else {
    inSCE <- temp.SCE
  }
  return(inSCE)
}




#' @title Subset a SingleCellExperiment object by columns
#' @description Used to peform subsetting of a
#' \linkS4class{SingleCellExperiment} object using a variety of methods that
#' indicate the correct columns to keep. The various methods,
#' \code{index}, \code{bool}, and \code{colData}, can be used in conjunction
#' with one another.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param index Integer vector. Vector of indicies indicating which columns
#' to keep. If \code{NULL}, this will not be used for subsetting.
#' Default \code{NULL}.
#' @param bool Boolean vector. Vector of \code{TRUE} or \code{FALSE}
#' indicating which columns should be kept. Needs to be the same length as the
#' number of columns in \code{inSCE}. If \code{NULL}, this will not be used
#' for subsetting. Default \code{NULL}.
#' @param colData Character. An expression that will identify a subset of
#' columns using variables found in the \code{colData} of \code{inSCE}.
#' For example, if \code{x} is a numeric vector in \code{colData},
#' then \code{"x < 5"} will return all columns with x less than 5.
#' Single quotes should be used for character strings. For example,
#' \code{"y == 'yes'"} will return all columns where y is "yes".
#' Multiple expressions can be evaluated by placing them in a vector.
#' For example \code{c("x < 5", "y =='yes'")} will apply both operations for
#' subsetting. If \code{NULL}, this will not be used for subsetting.
#' Default \code{NULL}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object that has
#' been subsetted by colData.
#' @author Joshua D. Campbell
#' @examples
#' data(scExample)
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' @export
subsetSCECols <- function(inSCE, index = NULL, bool = NULL, colData = NULL) {

  if(is.null(index) & is.null(bool) & is.null(colData)) {
    stop("At least one of 'index', 'bool', or 'colData' must be supplied.")
  }
  final.ix <- rep(FALSE, ncol(inSCE))

  # Parse index containing integers
  if(!is.null(index)) {
    if(min(index) < 1 | max(index) > ncol(inSCE)) {
      stop("'index' must contain integers between 1 and the number of columns ",
           "in 'inSCE': ", ncol(inSCE))
    }
    final.ix[index] <- TRUE
  }

  # Parse Boolean vector
  if(!is.null(bool)) {
    if(length(bool) != ncol(inSCE) | !is.logical(bool)) {
      stop("'bool' must be a logical vector the same length as the number of ",
           "colmns in 'inSCE': ", ncol(inSCE))
    }
    final.ix[bool] <- TRUE
  }

  # Parse expressions for colData variables
  if(!is.null(colData)) {
    for(i in seq_along(colData)) {
      temp <- eval(parse(text = as.character(colData[i])),
                   envir = as.data.frame(colData(inSCE)))
      if(length(temp) != ncol(inSCE) | !is.logical(temp)) {
        stop("The expression ", colData[i], " did not produce a boolean ",
             "vector the same length as the number of columns in 'inSCE'. ",
             "Please ensure that the spelling of the variable you are ",
             "trying to use matches one of the column names in ",
             "'colData(inSCE)' and that your expression is valid.")
      }
      final.ix[temp] = TRUE
    }
  }

  inSCE <- inSCE[,final.ix]
  return(inSCE)
  }
