#remove tag from metadata if NULL is set to assay

#' expDeleteDataTag
#' Remove tag against an input data from the stored tag information in the metadata of the input object.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assay Name of the assay or the data item against which a tag should be removed.
#' @return The input \code{SingleCellExperiment} object with tag information removed from the metadata slot.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- expSetDataTag(sce, "raw", "counts")
#' sce <- expDeleteDataTag(sce, "counts")
expDeleteDataTag <- function(inSCE, assay){
  if(!is.null(S4Vectors::metadata(inSCE)$assayType)){
    tbl <- S4Vectors::metadata(inSCE)$assayType
    tbl <- tbl %>% dplyr::filter(!.data$assayName %in% assay)
    S4Vectors::metadata(inSCE)$assayType <- tbl
  }
  return(inSCE)
}

#' expSetDataTag
#' Set tag to an assay or a data item in the input SCE object.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayType Specify a \code{character(1)} value as a tag that should be set against a data item.
#' @param assays Specify name(s) \code{character()} of data item(s) against which the tag should be set.
#' @return The input \code{SingleCellExperiment} object with tag information stored in the metadata slot.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- expSetDataTag(sce, "raw", "counts")
expSetDataTag <- function(inSCE, assayType, assays){
  tbl <- NULL
  if(is.null(S4Vectors::metadata(inSCE)$assayType)){
    tbl <- tibble::tibble(assayTag = assayType, assayName = assays)
  }
  else{
    tbl <- S4Vectors::metadata(inSCE)$assayType
    tbl <- rbind(tbl, tibble::tibble(assayTag = assayType, assayName = assays))
  }
  S4Vectors::metadata(inSCE)$assayType <- tbl
  return(inSCE)
}

#' expTaggedData
#' Returns a list of names of data items from the 
#' input \code{SingleCellExperiment} object based upon the input parameters.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param tags A \code{character()} value indicating if the data items should 
#' be returned separated by the specified tags. Default is \code{NULL} 
#' indicating that returned names of the data items are simply returned as a 
#' list with default tag as "uncategorized".
#' @param redDims A \code{logical} value indicating if \code{reducedDims} 
#' should be returned as well separated with 'redDims' tag.
#' @param recommended A \code{character()} vector indicating the tags that 
#' should be displayed as recommended. Default is \code{NULL}.
#' @param showTags A \code{logical} value indicating if the tags should be 
#' shown. If \code{FALSE}, output is just a simple list, not separated by tags.
#' @return A \code{list} of names of data items specified by the other 
#' parameters.
#' @importFrom stats filter
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- expSetDataTag(sce, "raw", "counts")
#' tags <- expTaggedData(sce)
expTaggedData <- function(inSCE, tags = NULL, redDims = FALSE, recommended = NULL, showTags = TRUE){
  namedList <- NULL
  tbl <- S4Vectors::metadata(inSCE)$assayType 
  
  if(!is.null(tags)){
    tbl <- tbl %>% dplyr::filter(.data$assayTag %in% tags)
  }
  
  if(redDims){
    tbl <- rbind(tbl, tibble::tibble(assayTag = "redDims", assayName = reducedDimNames(inSCE)))
  }
  
  if(!is.null(recommended)){
    recIx <- which(tbl$assayTag %in% recommended)
    if(length(recIx) > 0){
      tbl[recIx, ]$assayTag <- paste0(tbl$assayTag[recIx], " (recommended)")
      tbl <- rbind(tbl[recIx, ], tbl[-recIx, ])
      tbl$assayTag <- factor(tbl$assayTag, levels=unique(tbl$assayTag))
    }
  }
  
  if(!showTags){
    namedList <- as.character(tbl$assayName)
  }
  else{
    namedList <- with(tbl, split(tbl$assayName, tbl$assayTag))
    namedList <- lapply(namedList, vapply, list, list(length(namedList)))
  }

  return(namedList)
}

setClassUnion("CharacterOrNullOrMissing", c("character", "NULL", "missing"))
#' expData
#' Store data items using tags to identify the type of data item stored. To be used as a replacement for assay<- setter function but with additional parameter to set a tag to a data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @param tag Specify the tag to store against the input assay. Default is \code{NULL}, which will set the tag to "uncategorized".
#' @param altExp A \code{logical} value indicating if the input assay is a \code{altExp} or a subset assay.
#' @param value An input matrix-like value to store in the SCE object.
#' @return A \code{SingleCellExperiment} object containing the newly stored data.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' mat <- expData(sce, "counts")
#' expData(sce, "counts", tag = "raw") <- mat
setGeneric(name = "expData<-",
           function(inSCE, assayName, tag = NULL, altExp = FALSE, value)
             SummarizedExperiment::`assay<-`(x = inSCE,
                                             i = assayName,
                                             value = value)
)

#' expData
#' Store data items using tags to identify the type of data item stored. To be used as a replacement for assay<- setter function but with additional parameter to set a tag to a data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @param tag Specify the tag to store against the input assay. Default is \code{NULL}, which will set the tag to "uncategorized".
#' @param altExp A \code{logical} value indicating if the input assay is a \code{altExp} or a subset assay.
#' @param value An input matrix-like value to store in the SCE object.
#' @return A \code{SingleCellExperiment} object containing the newly stored data.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' mat <- expData(sce, "counts")
#' expData(sce, "counts", tag = "raw") <- mat
setMethod(f = "expData<-",
          signature = signature(
            inSCE = "ANY",
            assayName = "character",
            tag = "CharacterOrNullOrMissing",
            altExp = "logical"),
          definition = function(inSCE, assayName, tag = NULL, altExp = FALSE, value){
            if(!is.null(value)){
              if(is.null(tag)
                 || missing(tag)){
                inSCE <- expSetDataTag(
                  inSCE = inSCE,
                  assayType = "uncategorized",
                  assays = assayName
                )
              }
              else{
                inSCE <- expSetDataTag(
                  inSCE = inSCE,
                  assayType = tag,
                  assays = assayName
                )
              }
            }
            else{
              inSCE <- expDeleteDataTag(
                inSCE = inSCE,
                assay = assayName
              )
            }
            if(altExp){
              if (inherits(value, "SummarizedExperiment")) {
                altExp(inSCE, assayName) <- value
              } else {
                altExp(inSCE, assayName) <- SingleCellExperiment(list(counts = value))
                SummarizedExperiment::assayNames(altExp(inSCE, assayName)) <- assayName
              }
            }
            else{
                inSCE <- methods::callNextMethod()
            }
            return(inSCE)
          }
)

#' expData
#' Get data item from an input \code{SingleCellExperiment} object. The data item can be an \code{assay}, \code{altExp} (subset) or a \code{reducedDim}, which is retrieved based on the name of the data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the data item to retrieve.
#' @return Specified data item.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' mat <- expData(sce, "counts")
setGeneric(name = "expData",
           function(inSCE, assayName)
             SummarizedExperiment::assay(x = inSCE,
                                         i = assayName)
)

#' expData
#' Get data item from an input \code{SingleCellExperiment} object. The data item can be an \code{assay}, \code{altExp} (subset) or a \code{reducedDim}, which is retrieved based on the name of the data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the data item to retrieve.
#' @return Specified data item.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' mat <- expData(sce, "counts")
setMethod(f = "expData",
          signature = signature(inSCE = "ANY", assayName = "character"),
          definition = function(inSCE,  assayName){
            result <- NULL
            if(assayName %in% altExpNames(inSCE)){
              result <- altExp(inSCE, assayName)
              if(nrow(inSCE)<=nrow(result)
                 && ncol(inSCE) <= ncol(result)){
                inSCE <- result[rownames(inSCE), colnames(inSCE)]
              }
              else{
                inSCE <- result
              }
              result <- methods::callNextMethod()
            }
            else if (assayName %in% reducedDimNames(inSCE)){
              result <- reducedDim(inSCE, assayName)
            }
            else{
              result <- methods::callNextMethod()
            }
           return(result)
          }
)

#' expDataNames
#' Get names of all the data items in the input \code{SingleCellExperiment} object including assays, altExps and reducedDims.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @return A combined \code{vector} of \code{assayNames}, \code{altExpNames} and \code{reducedDimNames}.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' expDataNames(sce)
setGeneric(name = "expDataNames",
           function(inSCE)
             SummarizedExperiment::assayNames(x = inSCE)
)

#' expDataNames
#' Get names of all the data items in the input \code{SingleCellExperiment} object including assays, altExps and reducedDims.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @return A combined \code{vector} of \code{assayNames}, \code{altExpNames} and \code{reducedDimNames}.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' expDataNames(sce)
setMethod(f = "expDataNames",
          signature = signature(inSCE = "ANY"),
          definition = function(inSCE){
            result <- c(methods::callNextMethod(), altExpNames(inSCE), reducedDimNames(inSCE))
            return(result)
          }
)

# This utility function removes the redDims from untagged assays for use in expTaggedData function.
.filterRedDims <- function(inSCE, untaggedAssays){
  notRedDims <- untaggedAssays[!untaggedAssays %in% reducedDimNames(inSCE)]
  return(notRedDims)
}
