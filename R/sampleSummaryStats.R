#' @export
setGeneric("sampleSummaryStatsTable", function(inSCE, ...) standardGeneric("sampleSummaryStatsTable"))

#' @export
setGeneric("sampleSummaryStatsTable<-", function(inSCE, ..., value) standardGeneric("sampleSummaryStatsTable<-"))

#' @export
setMethod("sampleSummaryStatsTable", "SingleCellExperiment", function(inSCE, statsName, ...){
    return(inSCE@metadata$sctk$sampleSummary[[statsName]])
})

#' @rdname sampleSummaryStatsTable
#' @title Stores and returns table of SCTK QC outputs to metadata.
#' @description  Stores and returns table of QC metrics generated from
#'  QC algorithms within the metadata slot of the SingleCellExperiment object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' \link{assay} data and/or \link{colData} data. Required.
#' @param statsName A \code{character} value indicating the slot
#' that stores the stats table within the metadata of the
#' SingleCellExperiment object. Required.
#' @return A matrix/array object. Contains a summary table for QC statistics
#' generated from SingleCellTK.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- sampleSummaryStats(sce, simple = TRUE)
#' sampleSummaryStatsTable(sce, statsName = "sctk_qc")
#' @export
setReplaceMethod("sampleSummaryStatsTable", c("SingleCellExperiment", "ANY"), function(inSCE, statsName, ..., value) {
    inSCE@metadata$sctk$sampleSummary[[statsName]] <- value
    return(inSCE)
})

.sampleSummaryStats <- function(inSCE, colName = "Total",
                                simple = TRUE){

    metrics <- c("Number of Cells")
    values <- c(as.integer(ncol(inSCE)))

    if ("sum" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(metrics, "Mean counts", "Median counts")
        values <- c(values, mean(inSCE$sum),
                    stats::median(inSCE$sum))
    }

    if ("detected" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(
            metrics, "Mean features detected",
            "Median features detected"
        )
        values <- c(values, mean(inSCE$detected),
                    stats::median(inSCE$detected))
    }

    if(simple != TRUE){
        if ("dropletUtils_barcodeRank_knee" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(
                metrics, "BarcodeRank - Number of libraries above knee point"
            )
            values <- c(values, sum(SummarizedExperiment::colData(inSCE)$dropletUtils_BarcodeRank_Knee))
        }

        if ("dropletUtils_barcodeRank_inflection" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(
                metrics, "BarcodeRank - Number of libraries above inflection point"
            )
            values <- c(values, sum(SummarizedExperiment::colData(inSCE)$dropletUtils_BarcodeRank_Inflection))
        }

        if ("scrublet_call" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(
                metrics, "Scrublet - Number of doublets",
                "Scrublet - Percentage of doublets"
            )
            values <- c(values, sum(inSCE$scrublet_call == TRUE),
                        signif(sum(inSCE$scrublet_call == TRUE) / length(inSCE$scrublet_call) * 100, 3)
            )
        }

        if ("scDblFinder_doublet_call" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(metrics, "scDblFinder - Number of doublets",
                         "scDblFinder - Percentage of doublets")
            values <- c(values, sum(inSCE$scDblFinder_doublet_call == "Doublet"),
                        signif(sum(inSCE$scDblFinder_doublet_call == "Doublet")/length(inSCE$scDblFinder_doublet_call) * 100, 3))
        }

        if (any(grepl("doubletFinder_doublet_label_resolution",
                      colnames(SummarizedExperiment::colData(inSCE))))) {
            dfIx <- grep("doubletFinder_doublet_label_resolution",
                         colnames(SummarizedExperiment::colData(inSCE)))
            for (ix in dfIx) {
                metrics <- c(metrics,
                             paste("DoubletFinder - Number of doublets, Resolution",
                                   gsub("doubletFinder_doublet_label_resolution_",
                                        "",
                                        colnames(SummarizedExperiment::colData(inSCE))[ix])),
                             paste("DoubletFinder - Percentage of doublets, Resolution",
                                   gsub("doubletFinder_doublet_label_resolution_",
                                        "",
                                        colnames(SummarizedExperiment::colData(inSCE))[ix])))
                values <- c(values,sum(SummarizedExperiment::colData(inSCE)[, ix] == "Doublet"),
                            signif(sum(SummarizedExperiment::colData(inSCE)[, ix] == "Doublet") / length(SummarizedExperiment::colData(inSCE)[, ix]) * 100, 3))
            }
        }

        if("scds_cxds_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "CXDS - Number of doublets",
                         "CXDS - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_cxds_call == TRUE),
                        signif(sum(inSCE$scds_cxds_call == TRUE)/length(inSCE$scds_cxds_call) * 100, 3))
        }

        if("scds_bcds_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "BCDS - Number of doublets",
                         "BCDS - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_bcds_call == TRUE),
                        signif(sum(inSCE$scds_bcds_call == TRUE)/length(inSCE$scds_bcds_call) * 100, 3))
        }

        if("scds_hybrid_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "SCDS Hybrid - Number of doublets",
                         "SCDS Hybrid - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_hybrid_call == TRUE),
                        signif(sum(inSCE$scds_hybrid_call == TRUE)/length(inSCE$scds_hybrid_call) * 100, 3))
        }

        if("decontX_clusters" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "DecontX - Mean contamination",
                         "DecontX - Median contamination")
            values <- c(values, signif(mean(inSCE$decontX_contamination), 3),
                        signif(stats::median(inSCE$decontX_contamination), 3))
        }
    }

    df <- matrix(values)
    rownames(df) <- metrics
    colnames(df) <- colName
    return(df)
}

#' @title Generate table of SCTK QC outputs.
#' @description  Creates a table of QC metrics generated from
#'  QC algorithms, which is stored within the metadata slot of the
#'  input SingleCellExperiment object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' \link{assay} data and/or \link{colData} data. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#'  'counts'.
#' @param simple Boolean. Indicates whether to generate a table of only
#' basic QC stats (ex. library size), or to generate a summary table of all
#' QC stats stored in the inSCE.
#' @return A SingleCellExperiment object with a summary table for QC statistics
#' in the `sampleSummary` slot of metadata.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- sampleSummaryStats(sce, simple = TRUE)
#' sampleSummaryStatsTable(sce)
#' @importFrom magrittr %>%
#' @export
sampleSummaryStats <- function(inSCE,
                               sample = NULL,
                               useAssay = "counts",
                               simple = TRUE){

    if (!is.null(sample)) {
        if(length(sample) == 1){
            if(length(which(colnames(
                SummarizedExperiment::colData(inSCE)) == sample))){
                sample = colnames(SummarizedExperiment::colData(inSCE))
            }else{
                stop("'sample' must be either stored in the colData",
                     "or the same length as the number of columns in 'inSCE'")
            }

        }else if (length(sample) != ncol(inSCE)) {
            stop(
                "'sample' must be the same length as the number",
                " of columns in 'inSCE'"
            )
        }

    } else {
        sample <- rep("Sample", ncol(inSCE))
    }

    samples <- unique(sample)

    if(any(!c("sum", "detected") %in% colnames(SummarizedExperiment::colData(inSCE)))){
        inSCE <- scater::addPerCellQC(inSCE)
    }

    # if(simple == FALSE){
    dfTableAll <- .sampleSummaryStats(inSCE = inSCE, simple = simple)

    if(length(samples) > 1){
        forTable <- lapply(samples, function(x) {
            sampleInd <- which(sample == x)
            sampleSub <- sample[sampleInd]
            inSCESub <- inSCE[, sampleInd]
            df <- .sampleSummaryStats(inSCE = inSCESub, colName = x,
                                      simple = simple)
            return(df)
        })
        dfTableSample <- do.call(cbind, forTable)
        dfTableRes <- cbind(dfTableSample, dfTableAll)
    }else{
        dfTableRes <- dfTableAll
    }

    dfTableRes <- as.data.frame(dfTableRes)

    dfTableRes <- apply(dfTableRes, seq(2), function(x){
        return((signif(x,5)))
    })

    sampleSummaryStatsTable(inSCE, statsName = "sctk_qc") <- dfTableRes
    return(inSCE)
}