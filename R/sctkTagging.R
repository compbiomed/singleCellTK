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

#' @export
sctkSetTag <- function(inSCE, assayType, assays){
  if(!assays %in% S4Vectors::metadata(inSCE)$assayType[[assayType]]){
    S4Vectors::metadata(inSCE)$assayType[[assayType]] <- base::append(S4Vectors::metadata(inSCE)$assayType[[assayType]], assays)
  }
  #S4Vectors::metadata(inSCE)$assayType[[assayType]] <- assays
  return(inSCE)
}

#' @export
sctkSetTagExternal <- function(inSCE, assayType, assays){
  # S4Vectors::metadata(inSCE)$assayType[[assayType]] <- base::append(S4Vectors::metadata(inSCE)$assayType[[assayType]], assays)
  S4Vectors::metadata(inSCE)$assayType[[assayType]] <- assays
  return(inSCE)
}

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
#' @param value Input matrix-type assay to store.
#' @export
setGeneric(name = "sctkAssay<-", 
           function(inSCE, assayName, tag = NULL, value) 
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
#' @param value Input matrix-type assay to store.
#' @export
setMethod(f = "sctkAssay<-", 
          signature = signature(inSCE = "ANY", assayName = "character", tag = "CharacterOrNullOrMissing"),
          definition = function(inSCE,  assayName, tag = NULL, value){
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
            methods::callNextMethod()
          }
)

#' sctkAltExp
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param tag Specify the tag to store against the input assay. Default
#'  is \code{NULL}, which will set the tag to 'uncategorized'.
#' @param value Input matrix-type assay to store.
#' @export
setGeneric(name = "sctkAltExp<-", 
           function(inSCE, e, tag = NULL, value) 
             SingleCellExperiment::`altExp<-`(x = inSCE, 
                                              e = e,
                                             value = value)
)

#' sctkAltExp
#' Store assays using tags to identify the type of assay stored. To be used
#' within the singleCellTK as a replacement for assay<- setter function.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param e altExp Name
#' @param tag Specify the tag to store against the input assay. Default
#'  is \code{NULL}, which will set the tag to 'uncategorized'.
#' @param value Input matrix-type assay to store.
#' @export
setMethod(f = "sctkAltExp<-", 
          signature = signature(inSCE = "ANY", e = "character", tag = "CharacterOrNullOrMissing"),
          definition = function(inSCE, e, tag = NULL, value){
            if(!is.null(value)){
              if(is.null(tag)
                 || missing(tag)){
                inSCE <- sctkSetTag(
                  inSCE = inSCE, 
                  assayType = "altExp", 
                  assays = e
                )
              }
              else{
                inSCE <- sctkSetTag(
                  inSCE = inSCE, 
                  assayType = tag, 
                  assays = e
                )
              }
            }
            else{
              inSCE <- sctkDeleteTag(
                inSCE = inSCE,
                assay = e
              )
            }
            methods::callNextMethod()
          }
)
