#remove tag from metadata if NULL is set to assay
#' @export
sctkDeleteTag <- function(inSCE, assay){
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

#set a tag to an assay(s), used in sctkassay
#' @export
sctkSetTag <- function(inSCE, assayType, assays){
  if(!assays %in% S4Vectors::metadata(inSCE)$assayType[[assayType]]){
    S4Vectors::metadata(inSCE)$assayType[[assayType]] <- base::append(S4Vectors::metadata(inSCE)$assayType[[assayType]], assays)
  }
  #S4Vectors::metadata(inSCE)$assayType[[assayType]] <- assays
  return(inSCE)
}

#when you upload the file and it already has assays stored, then give a tag to all assays (in import code)
#' @export
sctkSetTagExternal <- function(inSCE, assayType, assays){
  # S4Vectors::metadata(inSCE)$assayType[[assayType]] <- base::append(S4Vectors::metadata(inSCE)$assayType[[assayType]], assays)
  S4Vectors::metadata(inSCE)$assayType[[assayType]] <- assays
  return(inSCE)
}

#returns all assays names with their tags (if no tags, then returned as uncategorized)
#' @export
getAssays <- function(inSCE){
  retList <- list()
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
  return(retList)
}

#returns assaynames against a selected tag
#' @export
sctkGetTag <- function(inSCE, assayType){
  retList <- list()
  if(!is.null(S4Vectors::metadata(inSCE)$assayType)){
    for(i in seq(assayType)){
      if(!is.null(S4Vectors::metadata(inSCE)$assayType[[assayType[i]]])){
        if(length(S4Vectors::metadata(inSCE)$assayType[[assayType[i]]]) == 1){
          #doing this because of how selectInput named list works, otherwise not needed
          retList[[assayType[i]]] <- list(S4Vectors::metadata(inSCE)$assayType[[assayType[i]]])
        }
        else{
          retList[[assayType[i]]] <- S4Vectors::metadata(inSCE)$assayType[[assayType[i]]]
        }
      }
    }
  }
  else{
    S4Vectors::metadata(inSCE)$assayType[["uncategorized"]] <- SummarizedExperiment::assayNames(inSCE)
    retList[["uncategorized"]] <- S4Vectors::metadata(inSCE)$assayType[["uncategorized"]]
  }
  return(retList)
}

setClassUnion("CharacterOrNullOrMissing", c("character", "NULL", "missing"))
#' sctkAssay
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @param tag Specify the tag to store against the input assay. Default
#'  is \code{NULL}, which will set the tag to 'uncategorized'.
#' @param altExp Logical value
#' @param value Input matrix-type assay to store.
#' @export
setGeneric(name = "sctkAssay<-", 
           function(inSCE, assayName, tag = NULL, altExp = FALSE, value) 
             SummarizedExperiment::`assay<-`(x = inSCE, 
                                             i = assayName, 
                                             value = value)
)

#' sctkAssay
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @param tag Specify the tag to store against the input assay. Default
#'  is \code{NULL}, which will set the tag to 'uncategorized'.
#' @param altExp Logical value
#' @param value Input matrix-type assay to store.
#' @export
setMethod(f = "sctkAssay<-", 
          signature = signature(inSCE = "ANY", assayName = "character", tag = "CharacterOrNullOrMissing", altExp = "logical"),
          definition = function(inSCE,  assayName, tag = NULL, altExp = FALSE, value){
            if(!is.null(value)){
              if(is.null(tag)
                 || missing(tag)){
                inSCE <- sctkSetTag(
                  inSCE = inSCE, 
                  assayType = "uncategorized", 
                  assays = assayName
                )
              }
              else{
                inSCE <- sctkSetTag(
                  inSCE = inSCE, 
                  assayType = tag, 
                  assays = assayName
                )
              }
            }
            else{
              inSCE <- sctkDeleteTag(
                inSCE = inSCE,
                assay = assayName
              )
            }
            if(altExp){
              altExp(inSCE, assayName) <- SingleCellExperiment(list(counts = value))
              assayNames(altExp(inSCE, assayName)) <- assayName
            }
            else{
              inSCE <- methods::callNextMethod()
            }
            return(inSCE)
          }
)

#' sctkAssay
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @export
setGeneric(name = "sctkAssay", 
           function(inSCE, assayName) 
             SummarizedExperiment::assay(x = inSCE,
                                         i = assayName)
)

#' sctkAssay
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param assayName Specify the name of the input assay.
#' @export
setMethod(f = "sctkAssay", 
          signature = signature(inSCE = "ANY", assayName = "character"),
          definition = function(inSCE,  assayName){
            result <- NULL
            if(assayName %in% altExpNames(inSCE)){
              result <- altExp(inSCE, assayName)
              result <- assay(result, assayName)
            }
            else{
              result <- methods::callNextMethod()
            }
           return(result) 
          }
)

#' sctkAssayNames
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @export
setGeneric(name = "sctkAssayNames", 
           function(inSCE) 
             SummarizedExperiment::assayNames(x = inSCE)
)

#' sctkAssayNames
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @export
setMethod(f = "sctkAssayNames", 
          signature = signature(inSCE = "ANY"),
          definition = function(inSCE){
            result <- c(altExpNames(inSCE), methods::callNextMethod())
            return(result)
          }
)
