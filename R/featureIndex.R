#' @title Retrieve row index for a set of features
#' @description This will return indices of features among the rownames
#' or rowData of a data.frame, matrix, or a \linkS4class{SummarizedExperiment}
#' object including a \linkS4class{SingleCellExperiment}.
#' Partial matching (i.e. grepping) can be used by setting
#' \code{exactMatch = FALSE}.
#' @param features Character vector of feature names to find in the rows of
#' \code{inSCE}.
#' @param inSCE A data.frame, matrix, or \linkS4class{SingleCellExperiment}
#' object to search.
#' @param by Character. Where to search for features in \code{inSCE}. If set to
#' \code{"rownames"} then the features will be searched for among
#' \code{rownames(inSCE)}. If \code{inSCE} inherits from class
#' \linkS4class{SummarizedExperiment}, then \code{by} can be one of the
#' fields in the row annotation data.frame (i.e. one of
#' \code{colnames(rowData(inSCE))}).
#' @param exactMatch Boolean. Whether to only identify exact matches
#' or to identify partial matches using \code{\link{grep}}.
#' @param removeNA Boolean. If set to \code{FALSE}, features not found in
#' \code{inSCE} will be given \code{NA} and the returned vector will be the same
#' length as \code{features}. If set to \code{TRUE}, then the \code{NA}
#' values will be removed from the returned vector. Default \code{FALSE}.
#' @param errorOnNoMatch Boolean. If \code{TRUE}, an error will be given if
#' no matches are found. If \code{FALSE}, an empty vector will be returned if 
#' \code{removeNA} is set to \code{TRUE} or a vector of \code{NA} if 
#' \code{removeNA} is set to \code{FALSE}. Default \code{TRUE}.
#' @param warningOnPartialMatch Boolean. If \code{TRUE}, a warning will be
#' given if some of the entries in \code{features} were not found in 
#' \code{inSCE}. The warning will list the features not found.
#' Default \code{TRUE}.
#' @return A vector of row indices for the matching features in \code{inSCE}.
#' @author Yusuke Koga, Joshua D. Campbell
#' @seealso '\link[scater]{retrieveFeatureInfo}' from package \code{'scater'}
#' and \code{link{regex}} for how to use regular expressions when
#' \code{exactMatch = FALSE}.
#' @examples
#' data(scExample)
#' ix <- featureIndex(features = c("MT-CYB", "MT-ND2"),
#'                              inSCE = sce,
#'                              by = "feature_name")
#' @export
featureIndex <- function(features, inSCE,
                                 by = "rownames",
                                 exactMatch = TRUE,
                                 removeNA = FALSE,
                                 errorOnNoMatch = TRUE,
                                 warningOnPartialMatch = TRUE) {
  
  # Extract vector to search through
  if (by == "rownames") {
    if (is.null(rownames(inSCE))) {
      stop("'rownames' of 'inSCE' are 'NULL'. Please set 'rownames' or change",
           " 'by' to search a different column in 'inSCE'.")
    }
    search <- rownames(inSCE)
  } else if (length(ncol(inSCE)) > 0) {
    if (inherits(inSCE, "SummarizedExperiment")) {
      if (!(by %in% colnames(rowData(inSCE)))) {
        stop("'by' is not a column in 'rowData(inSCE)'.")
      }
      search <- rowData(inSCE)[, by]
    } else {
      if (!(by %in% colnames(inSCE))) {
        stop("'by' is not a column in 'inSCE'.")
      }
      search <- inSCE[, by]
    }
  } else {
    search <- as.character(inSCE)
  }
  
  # Match each element of 'pattern' in vector 'search'
  if (!isTRUE(exactMatch)) {
    featuresIndices <- rep(NA, length(features))
    for (i in seq_along(features)) {
      g <- grep(features[i], search)
      if (length(g) == 1) {
        featuresIndices[i] <- g
      } else if (length(g) > 1) {
        warning(
          "Feature '", features[i], "' matched multiple items in '",
          by, "': ", paste(search[g], collapse = ","),
          ". Only the first match will be selected."
        )
        featuresIndices[i] <- g[1]
      }
    }
  } else {
    featuresIndices <- match(features, search)
  }
  
  if (sum(is.na(featuresIndices)) > 0) {
    if (sum(is.na(featuresIndices)) == length(features)
        & isTRUE(errorOnNoMatch)) {
      if (isTRUE(exactMatch)) {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'inSCE'. ",
          "Check the spelling or try setting",
          " 'exactMatch = FALSE'."
        )
      } else {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'inSCE'. ",
          "Check the spelling and make sure 'by' is set",
          " to the appropriate place in 'inSCE'."
        )
      }
    }
    if(isTRUE(warningOnPartialMatch)){
      warning(
        "The following features were not present in 'inSCE': ",
        paste(features[which(is.na(featuresIndices))],
              collapse = ","
        )
      )
    }
  }
  
  if (isTRUE(removeNA)) {
    featuresIndices <- featuresIndices[!is.na(featuresIndices)]
  }
  return(featuresIndices)
}
