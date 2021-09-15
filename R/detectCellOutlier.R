#' @title Detecting outliers within the SingleCellExperiment object.
#' @description A wrapper function for \link[scater]{isOutlier}. Identify
#'  outliers from numeric vectors stored in the SingleCellExperiment object.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample A single character specifying a name that can be found in
#' \code{colData(inSCE)} to directly use the cell annotation; or a character
#' vector with as many elements as cells to indicates which sample each cell
#' belongs to. Default NULL. \link[celda]{decontX} will be run on cells from
#' each sample separately.
#' @param slotName Desired slot of SingleCellExperiment used for plotting. Possible
#'  options: "assays", "colData", "metadata", "reducedDims". Required.
#' @param itemName Desired vector within the slot used for plotting. Required.
#' @param nmads Integer. Number of median absolute deviation. Parameter may be
#'  adjusted for more lenient or stringent outlier cutoff. Default 3.
#' @param type Character. Type/direction of outlier detection; whether
#'  the lower/higher outliers should be detected, or both.
#'  Options are "both", "lower", "higher".
#' @param overwrite Boolean. If TRUE, and this function has previously generated
#'  an outlier decision on the same itemName, the outlier decision will be overwritten.
#'  Default TRUE.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  '' added to the
#'  \link{colData} slot. Additionally, the
#' decontaminated counts will be added as an assay called 'decontXCounts'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDecontX(sce[,sample(ncol(sce),20)])
#' sce <- detectCellOutlier(sce, slotName = "colData", sample = sce$sample,
#'  nmads = 4, itemName = "decontX_contamination", type = "both")
#' @export
detectCellOutlier <- function(inSCE, slotName, itemName, sample = NULL, nmads = 3,
                              type = "both", overwrite = TRUE){
    if (!slotName %in% c("rowData", "colData", "assays", "metadata", "reducedDims")) {
        stop("'slotName' must be a slotName within the SingleCellExperiment object.",
             "Please run 'methods::slot' if you are unsure the",
             "specified slotName exists.")
    }
    if(!is.numeric(nmads)){
        nmads = 3
    }
    argsList <- mget(names(formals()),sys.frame(sys.nframe()))
    outlierDF <- S4Vectors::DataFrame(row.names = colnames(inSCE),
                                      outlier = logical(ncol(inSCE)),
                                      threshold = numeric(ncol(inSCE)))
    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                 " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- inSCE[, sceSampleInd]

        value <- do.call(slotName, args = list(sceSample))[[itemName]]

        scaterRes <- scater::isOutlier(metric = value, type = type, nmads = nmads)
        outlierDF[sceSampleInd,1] <- scaterRes
        outlierDF[sceSampleInd,2] <- attr(scaterRes, "thresholds")[2]
    }

    colnames(outlierDF) <- paste(itemName, colnames(outlierDF), paste0("mad", nmads), type, sep = "_")
    if(overwrite){
        if(any(colnames(outlierDF) %in% colnames(colData(inSCE)))){
            replaceIx <- which(colnames(colData(inSCE)) %in% colnames(outlierDF))
            colData(inSCE)[,replaceIx] <- outlierDF
        }else{
            colData(inSCE) <- cbind(colData(inSCE), outlierDF)
        }
    }else{
        colData(inSCE) <- cbind(colData(inSCE), outlierDF)
    }
    inSCE@metadata$outlierDetection[[itemName]][[paste0("mad",nmads)]] <- argsList[-1]
    return(inSCE)
}