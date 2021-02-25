#remove tag from metadata if NULL is set to assay

#' expDeleteDataTag
#' Remove tag against an input data from the stored tag information in the metadata of the input object.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assay Name of the assay or the data item against which a tag should be removed.
#' @return The input \code{SingleCellExperiment} object with tag information removed from the metadata slot.
#' @export
#'
expDeleteDataTag <- function(inSCE, assay){
  for(i in seq(length(S4Vectors::metadata(inSCE)$assayType))){
    matchedIndex <- match(assay, S4Vectors::metadata(inSCE)$assayType[[i]])
    if(!is.na(matchedIndex)){
      if(length(S4Vectors::metadata(inSCE)$assayType[[i]]) == 1){
        S4Vectors::metadata(inSCE)$assayType[[i]] <- NULL
      }
      else{
        S4Vectors::metadata(inSCE)$assayType[[i]] <- S4Vectors::metadata(inSCE)$assayType[[i]][-matchedIndex]
      }
    }
  }
  return(inSCE)
}

#' expSetDataTag
#' Set tag to an assay or a data item in the input SCE object.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayType Specify a \code{character(1)} value as a tag that should be set against a data item.
#' @param assays Specify name(s) \code{character()} of data item(s) against which the tag should be set. 
#' @param append A \code{logical} value indicating if this assay should be appended to the object or overridden. Default value is \code{TRUE} indicating that it should be appended.
#' @return The input \code{SingleCellExperiment} object with tag information stored in the metadata slot.
#' @export
#'
expSetDataTag <- function(inSCE, assayType, assays, append = TRUE){
  if(append){
    if(!assays %in% S4Vectors::metadata(inSCE)$assayType[[assayType]]){
      S4Vectors::metadata(inSCE)$assayType[[assayType]] <- base::append(S4Vectors::metadata(inSCE)$assayType[[assayType]], assays)
    }
  }
  else{
    S4Vectors::metadata(inSCE)$assayType[[assayType]] <- assays
  }
  return(inSCE)
}

#' expTaggedData
#' Returns a list of names of data items from the input \code{SingleCellExperiment} object based upon the input parameters.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param tags A \code{character()} value indicating if the data items should be returned separated by the specified tags. Default is \code{NULL} indicating that returned names of the data items are simply returned as a list with default tag as "uncategorized".
#' @param redDims A \code{logical} value indicating if \code{reducedDims} should be returned as well separated with 'redDims' tag.
#' @return A \code{list} of names of data items specfied by the other parameters.
#' @export
#'
expTaggedData <- function(inSCE, tags = NULL, redDims = FALSE){
  retList <- list()
  if(is.null(tags)){
    if(!is.null(S4Vectors::metadata(inSCE)$assayType)){
      for(i in seq(S4Vectors::metadata(inSCE)$assayType)){
        if(!is.null(S4Vectors::metadata(inSCE)$assayType[[i]])){
          if(length(S4Vectors::metadata(inSCE)$assayType[[i]]) == 1){
            #doing this because of how selectInput named list works, otherwise not needed
            retList[[names(S4Vectors::metadata(inSCE)$assayType)[i]]] <- list(S4Vectors::metadata(inSCE)$assayType[[i]])
          }
          else{
            retList[[names(S4Vectors::metadata(inSCE)$assayType)[i]]] <- S4Vectors::metadata(inSCE)$assayType[[i]]
          }
        }
      }
    }
    else{
      S4Vectors::metadata(inSCE)$assayType[["uncategorized"]] <- SummarizedExperiment::assayNames(inSCE)
      retList[["uncategorized"]] <- S4Vectors::metadata(inSCE)$assayType[["uncategorized"]]
    }
  }
  else{
    if(!is.null(S4Vectors::metadata(inSCE)$assayType)){
      for(i in seq(tags)){
        if(!is.null(S4Vectors::metadata(inSCE)$assayType[[tags[i]]])){
          if(length(S4Vectors::metadata(inSCE)$assayType[[tags[i]]]) == 1){
            #doing this because of how selectInput named list works, otherwise not needed
            retList[[tags[i]]] <- list(S4Vectors::metadata(inSCE)$assayType[[tags[i]]])
          }
          else{
            retList[[tags[i]]] <- S4Vectors::metadata(inSCE)$assayType[[tags[i]]]
          }
        }
      }
    }
    else{
      S4Vectors::metadata(inSCE)$assayType[["uncategorized"]] <- SummarizedExperiment::assayNames(inSCE)
      retList[["uncategorized"]] <- S4Vectors::metadata(inSCE)$assayType[["uncategorized"]]
    }
  }
  if(redDims){
    if(length(reducedDimNames(inSCE)) > 0){
      if(length(reducedDimNames(inSCE)) == 1){
        retList[["redDims"]] <- list(reducedDimNames(inSCE))
      }
      else{
        retList[["redDims"]] <- reducedDimNames(inSCE)
      }
    }
  }
  return(retList)
}

setClassUnion("CharacterOrNullOrMissing", c("character", "NULL", "missing"))
#' expData
#' Store data items using tags to identify the type of data item stored. To be used as a replacement for assay<- setter function but with additional parameter to set a tag to a data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @param tag Specify the tag to store against the input assay. Default is \code{NULL}, which will set the tag to "uncategorized".
#' @param altExp A \code{logical} value indicating if the input assay is a \code{altExp} or a subset assay.
#' @param value An input matrix-like value to store in the SCE object.
#' @export
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
#' @export
setMethod(f = "expData<-", 
          signature = signature(
            inSCE = "ANY", 
            assayName = "character", 
            tag = "CharacterOrNullOrMissing", 
            altExp = "logical"),
          definition = function(inSCE,  assayName, tag = NULL, altExp = FALSE, value){
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
              altExp(inSCE, assayName) <- SingleCellExperiment(list(counts = value))
             SummarizedExperiment::assayNames(altExp(inSCE, assayName)) <- assayName
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
#' @export
setGeneric(name = "expData", 
           function(inSCE, assayName) 
             SummarizedExperiment::assay(x = inSCE,
                                         i = assayName)
)

#' expData
#' Get data item from an input \code{SingleCellExperiment} object. The data item can be an \code{assay}, \code{altExp} (subset) or a \code{reducedDim}, which is retrieved based on the name of the data item.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the data item to retrieve.
#' @export
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
setGeneric(name = "expDataNames", 
           function(inSCE) 
             SummarizedExperiment::assayNames(x = inSCE)
)

#' expDataNames
#' Get names of all the data items in the input \code{SingleCellExperiment} object including assays, altExps and reducedDims.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @return A combined \code{vector} of \code{assayNames}, \code{altExpNames} and \code{reducedDimNames}. 
#' @export
setMethod(f = "expDataNames", 
          signature = signature(inSCE = "ANY"),
          definition = function(inSCE){
            result <- c(altExpNames(inSCE), methods::callNextMethod())
            return(result)
          }
)
