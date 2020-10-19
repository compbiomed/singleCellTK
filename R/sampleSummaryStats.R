.sampleSummaryStats <- function(inSCE, colName = "Total",
                                simple = TRUE){

    metrics <- c("Number of Cells")
    values <- c(ncol(inSCE))

    if ("sum" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(metrics, "Mean counts", "Median counts")
        values <- c(values, signif(mean(inSCE$sum), 3),
                    signif(stats::median(inSCE$sum), 3))
    }

    if ("detected" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(
            metrics, "Mean features detected",
            "Median features detected"
        )
        values <- c(values, signif(mean(inSCE$detected), 3),
                    signif(stats::median(inSCE$detected), 3))
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

        if ("scran_doubletCells_score_log10" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(metrics, "DoubletCells - Doublet score outliers")
            values <- c(values, sum(scater::isOutlier(inSCE$scran_doubletCells_score_log10,
                                                  type = "higher"
            )))

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
#'  QC algorithms via either kable or csv file.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' \link[SummarizedExperiment]{assay} data and/or
#' \link[SummarizedExperiment]{colData} data. Required.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#'  'counts'.
#' @param simple Boolean. Indicates whether to generate a table of only
#' basic QC stats (ex. library size), or to generate a summary table of all
#' QC stats stored in the inSCE.
#' @return A matrix/array object.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sampleSummaryStats(sce, simple = TRUE)
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
    dfTableRes <- apply(dfTableRes, 1:2, function(x){
        if(grepl(as.character(x), "\\.0000")){
            return(as.integer(x))
        }else{
            return(as.numeric(x))
        }
    })

    return(dfTableRes)
}
