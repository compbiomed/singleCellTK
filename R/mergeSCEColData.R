#' @title Merging colData from two singleCellExperiment objects
#' @description Merges colData of the singleCellExperiment objects 
#'  obtained from the same dataset which contain differing colData. 
#'  (i.e. raw data and filtered data)  
#' @param sce1 SingleCellExperiment object. The function will output this
#'  singleCellExperiment object with a combined colData from sce1 and sce2.
#' @param sce2 SingleCellExperiment object. colData from this object
#'  will be merged with colData from sce1 and loaded into sce1.
#' @param id1 Character vector. Column in colData of sce1 that will be 
#'  used to combine sce1 and sce2. Default "column_name" 
#' @param id2 Character vector. Column in colData of sce2 that will be
#'  used to combine sce1 and sce2. Default "column_name"
#' @return SingleCellExperiment object containing combined colData from  
#'  both singleCellExperiment for samples in sce1.
#' @examples
#' sce <- mergeSCEColData(sce1, sce2, id1 = "column_name", id2 = "column_name")
#' @export
mergeSCEColData <- function(sce1, sce2, id1 = "column_name", id2 = "column_name"){
    #Subset coldata of sce2
    not.in.sce1 <- c(setdiff(names(SummarizedExperiment::colData(sce2)),
      names(SummarizedExperiment::colData(sce1))),id2)
    not.in.sce1 <- not.in.sce1[!is.null(not.in.sce1)]

    coldata.not.in.sce1 <- SummarizedExperiment::colData(sce2)[,c(not.in.sce1),
      drop = FALSE]

    coldata.sce1 <- SummarizedExperiment::colData(sce1)

    if(is.null(id1) | is.null(id2)){
        if(is.null(rownames(coldata.not.in.sce1))){
            stop("Unable to match between singleCellExperiment objects.
              Please define id1/id2 within the function,
              or assign a column name for the singleCellExperiment object.")
        }else{
            #Add placeholder for merging
            coldata.not.in.sce1$cell <- rownames(coldata.not.in.sce1)
            #Add placeholder for merging
            coldata.sce1$cell <- rownames(SummarizedExperiment::colData(sce1))
            id1 = "cell"
            id2 = "cell"
            placeholder = TRUE
        }
    }else{
        placeholder = FALSE
    }

    #Merge
    coldata.merge <- base::merge(coldata.sce1,
        coldata.not.in.sce1,
        all.x = TRUE,
        sort = FALSE,
        by.x = id1,
        by.y = id2)

    #Reorder to match counts matrix
    coldata.merge <- coldata.merge[match(colnames(SingleCellExperiment::counts(sce1)),
        coldata.merge[,id1]),]

    #Re-enter rownames back into coldata
    rownames(coldata.merge) <- coldata.merge[,id1]

    #Get rid of placeholder
    if(placeholder == TRUE){
    	  coldata.merge[,id1] <- NULL
    }

    SummarizedExperiment::colData(sce1) <- S4Vectors::DataFrame(coldata.merge)
    return(sce1)
}
