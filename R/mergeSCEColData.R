#' @title Merging colData from two singleCellExperiment objects
#' @description Merges colData of the singleCellExperiment objects 
#'  obtained from the same dataset which contain differing colData. 
#'  (i.e. raw data and filtered data)  
#' @param inSCE1 Input SingleCellExperiment object. The function will output this
#'  singleCellExperiment object with a combined colData from inSCE1 and inSCE2.
#' @param inSCE2 Input SingleCellExperiment object. colData from this object
#'  will be merged with colData from inSCE1 and loaded into inSCE1.
#' @param id1 Character vector. Column in colData of inSCE1 that will be 
#'  used to combine inSCE1 and inSCE2. Default "column_name" 
#' @param id2 Character vector. Column in colData of inSCE2 that will be
#'  used to combine inSCE1 and inSCE2. Default "column_name"
#' @return SingleCellExperiment object containing combined colData from  
#'  both singleCellExperiment for samples in inSCE1.
#' @examples
#' sce1 <- importCellRanger(
#'     cellRangerDirs = system.file("extdata/", package = "singleCellTK"),
#'     sampleDirs = "hgmm_1k_v3_20x20",
#'     sampleNames = "hgmm1kv3",
#'     dataType = "filtered")
#' data(scExample)
#' sce2 <- sce
#' sce <- mergeSCEColData(inSCE1 = sce1, inSCE2 = sce2, id1 = "column_name", id2 = "column_name")
#' @export
mergeSCEColData <- function(inSCE1, inSCE2, id1 = "column_name", id2 = "column_name"){
    not.in.sce1 <- c(setdiff(names(SummarizedExperiment::colData(inSCE2)),
      names(SummarizedExperiment::colData(inSCE1))),id2)
    not.in.sce1 <- not.in.sce1[!is.null(not.in.sce1)]

    coldata.not.in.sce1 <- SummarizedExperiment::colData(inSCE2)[,c(not.in.sce1),
      drop = FALSE]

    coldata.sce1 <- SummarizedExperiment::colData(inSCE1)

    if(is.null(id1) | is.null(id2)){
        if(is.null(rownames(coldata.not.in.sce1))){
            stop("Unable to match between singleCellExperiment objects.
              Please define id1/id2 within the function,
              or assign a column name for the singleCellExperiment object.")
        }else{
            coldata.not.in.sce1$cell <- rownames(coldata.not.in.sce1)
            coldata.sce1$cell <- rownames(SummarizedExperiment::colData(inSCE1))
            id1 <- "cell"
            id2 <- "cell"
            placeholder = TRUE
        }
    }else{
        placeholder <- FALSE
    }

    coldata.merge <- base::merge(coldata.sce1,
        coldata.not.in.sce1,
        all.x = TRUE,
        sort = FALSE,
        by.x = id1,
        by.y = id2)

    coldata.merge <- coldata.merge[match(colnames(SingleCellExperiment::counts(inSCE1)),
        coldata.merge[,id1]),]

    rownames(coldata.merge) <- coldata.merge[,id1]

    if(placeholder == TRUE){
    	  coldata.merge[,id1] <- NULL
    }

    SummarizedExperiment::colData(inSCE1) <- S4Vectors::DataFrame(coldata.merge)
    return(inSCE1)
}
