#' @rdname getSampleSummaryStatsTable
setMethod("getSampleSummaryStatsTable", "SingleCellExperiment", function(inSCE, statsName, ...){
    allStatsNames <- listSampleSummaryStatsTables(inSCE)
    if(!statsName %in% allStatsNames){
        stop(paste(statsName, "is not a table within the SingleCellExperiment object.",
                   "The following are the names of the tables stored:",
                   paste(allStatsNames, collapse = ",")))
    }else{
        return(inSCE@metadata$sctk$sample_summary[[statsName]])
    }
})

#' @rdname setSampleSummaryStatsTable<-
setReplaceMethod("setSampleSummaryStatsTable", c("SingleCellExperiment", "ANY"), function(inSCE, statsName, ..., value) {
    inSCE@metadata$sctk$sample_summary[[statsName]] <- value
    return(inSCE)
})

#' @rdname listSampleSummaryStatsTables
setMethod("listSampleSummaryStatsTables", "SingleCellExperiment", function(inSCE, ...){
    if(is.null(inSCE@metadata$sctk$sample_summary)){
        stop(paste("No sample-level QC tables are available.",
                     "Please try executing functions such as sampleSummaryStats first."))
    }else{
        allStatsNames <- names(inSCE@metadata$sctk$sample_summary)
        if(is.null(allStatsNames) || length(allStatsNames) == 0){
            stop(paste("No sample-level QC tables are available.",
                         "Please try executing functions such as sampleSummaryStats first."))
        }else{
            return(names(inSCE@metadata$sctk$sample_summary))
        }
    }
})

.sampleSummaryStats <- function(dataFrame, colName = "All Samples",
                                simple = TRUE){

    metrics <- c("Number of Cells")
    values <- c(as.integer(nrow(dataFrame)))

    if ("sum" %in% colnames(dataFrame)) {
        metrics <- c(metrics, "Mean counts", "Median counts")
        values <- c(values, mean(dataFrame$sum),
                    stats::median(dataFrame$sum))
    }

    if ("detected" %in% colnames(dataFrame)) {
        metrics <- c(
            metrics, "Mean features detected",
            "Median features detected"
        )
        values <- c(values, mean(dataFrame$detected),
                    stats::median(dataFrame$detected))
    }

    if(simple != TRUE){
        if ("dropletUtils_barcodeRank_knee" %in% colnames(dataFrame)) {
            metrics <- c(
                metrics, "BarcodeRank - Number of libraries above knee point"
            )
            values <- c(values, sum(dataFrame$dropletUtils_BarcodeRank_Knee))
        }

        if ("dropletUtils_barcodeRank_inflection" %in% colnames(dataFrame)) {
            metrics <- c(
                metrics, "BarcodeRank - Number of libraries above inflection point"
            )
            values <- c(values, sum(dataFrame$dropletUtils_BarcodeRank_Inflection))
        }

        if ("scrublet_call" %in% colnames(dataFrame)) {
            metrics <- c(
                metrics, "Scrublet - Number of doublets",
                "Scrublet - Percentage of doublets"
            )
            values <- c(values, sum(dataFrame$scrublet_call == TRUE),
                        signif(sum(dataFrame$scrublet_call == TRUE) / length(dataFrame$scrublet_call) * 100, 3)
            )
        }

        if ("scDblFinder_doublet_call" %in% colnames(dataFrame)) {
            metrics <- c(metrics, "scDblFinder - Number of doublets",
                         "scDblFinder - Percentage of doublets")
            values <- c(values, sum(dataFrame$scDblFinder_doublet_call == "Doublet"),
                        signif(sum(dataFrame$scDblFinder_doublet_call == "Doublet")/length(dataFrame$scDblFinder_doublet_call) * 100, 3))
        }

        if (any(grepl("doubletFinder_doublet_label_resolution",
                      colnames(dataFrame)))) {
            dfIx <- grep("doubletFinder_doublet_label_resolution",
                         colnames(dataFrame))
            for (ix in dfIx) {
                metrics <- c(metrics,
                             paste("DoubletFinder - Number of doublets, Resolution",
                                   gsub("doubletFinder_doublet_label_resolution_",
                                        "",
                                        colnames(dataFrame)[ix])),
                             paste("DoubletFinder - Percentage of doublets, Resolution",
                                   gsub("doubletFinder_doublet_label_resolution_",
                                        "",
                                        colnames(dataFrame)[ix])))
                values <- c(values,sum(dataFrame[, ix] == "Doublet"),
                            signif(sum(dataFrame[, ix] == "Doublet") / length(dataFrame[, ix]) * 100, 3))
            }
        }

        if("scds_cxds_call" %in% colnames(dataFrame)){
            metrics <- c(metrics, "CXDS - Number of doublets",
                         "CXDS - Percentage of doublets")
            values <- c(values, sum(dataFrame$scds_cxds_call == "Doublet"),
                        signif(sum(dataFrame$scds_cxds_call == "Doublet")/length(dataFrame$scds_cxds_call) * 100, 3))
        }

        if("scds_bcds_call" %in% colnames(dataFrame)){
            metrics <- c(metrics, "BCDS - Number of doublets",
                         "BCDS - Percentage of doublets")
            values <- c(values, sum(dataFrame$scds_bcds_call == "Doublet"),
                        signif(sum(dataFrame$scds_bcds_call == "Doublet")/length(dataFrame$scds_bcds_call) * 100, 3))
        }

        if("scds_hybrid_call" %in% colnames(dataFrame)){
            metrics <- c(metrics, "SCDS Hybrid - Number of doublets",
                         "SCDS Hybrid - Percentage of doublets")
            values <- c(values, sum(dataFrame$scds_hybrid_call == "Doublet"),
                        signif(sum(dataFrame$scds_hybrid_call == "Doublet")/length(dataFrame$scds_hybrid_call) * 100, 3))
        }

        if("decontX_clusters" %in% colnames(dataFrame)){
            metrics <- c(metrics, "DecontX - Mean contamination",
                         "DecontX - Median contamination")
            values <- c(values, signif(mean(dataFrame$decontX_contamination), 3),
                        signif(stats::median(dataFrame$decontX_contamination), 3))
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
#' @param statsName Character. The name of the slot that will store the
#' QC stat table. Default "qc_table".
#' @return A SingleCellExperiment object with a summary table for QC statistics
#' in the `sample_summary` slot of metadata.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- sampleSummaryStats(sce, simple = TRUE)
#' getSampleSummaryStatsTable(sce, statsName = "qc_table")
#' @importFrom magrittr %>%
#' @export
sampleSummaryStats <- function(inSCE,
                               sample = NULL,
                               useAssay = "counts",
                               simple = TRUE,
                               statsName = "qc_table"){

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
    qcDataFrame <- SummarizedExperiment::colData(inSCE)
    dfTableAll <- .sampleSummaryStats(dataFrame = qcDataFrame, simple = simple)

    if(length(samples) > 1){
        forTable <- lapply(samples, function(x) {
            sampleInd <- which(sample == x)
            sampleSub <- sample[sampleInd]
            inSCESub <- inSCE[, sampleInd]
            qcDataFrame <- SummarizedExperiment::colData(inSCESub)
            df <- .sampleSummaryStats(dataFrame = qcDataFrame, colName = x,
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

    setSampleSummaryStatsTable(inSCE, statsName = statsName) <- dfTableRes
    return(inSCE)
}
