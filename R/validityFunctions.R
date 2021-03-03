.checkSCEValidity <- function(inSCE){
  if(is.null(rownames(inSCE)))
    stop("Rownames of the input SCE object cannot be NULL.")
  if(is.null(colnames(inSCE)))
    stop("Colnames of the input SCE object cannot be NULL.")
}