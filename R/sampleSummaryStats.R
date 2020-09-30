.sampleSummaryStats <- function(inSCE, colName = "Total",
                                simple = TRUE){

    metrics <- c("Number of Cells")
    values <- c(ncol(inSCE))

    if ("sum" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(metrics, "Mean counts", "Median counts")
        values <- c(values, mean(inSCE$sum), stats::median(inSCE$sum))
    }

    if ("detected" %in% colnames(SummarizedExperiment::colData(inSCE))) {
        metrics <- c(
            metrics, "Mean features detected",
            "Median features detected"
        )
        values <- c(values, mean(inSCE$detected), stats::median(inSCE$detected))
    }

    if(simple != TRUE){
        if ("scrublet_call" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(
                metrics, "Scrublet - Number of doublets",
                "Scrublet - Percentage of doublets"
            )
            values <- c(
                values, sum(inSCE$scrublet_call == TRUE),
                sum(inSCE$scrublet_call == TRUE) / length(inSCE$scrublet_call) * 100
            )
        }

        if ("scran_doubletCells_score_log10" %in% colnames(SummarizedExperiment::colData(inSCE))) {
            metrics <- c(metrics, "DoubletCells - Doublet score outliers")
            values <- c(values, scater::isOutlier(inSCE$scran_doubletCells_logscore,
                                                  type = "higher"
            ))

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
                            sum(SummarizedExperiment::colData(inSCE)[, ix] == "Doublet") / length(SummarizedExperiment::colData(inSCE)[, ix]) * 100)
            }
        }

        if("scds_cxds_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "CXDS - Number of doublets",
                         "CXDS - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_cxds_call == TRUE),
                        sum(inSCE$scds_cxds_call == TRUE)/length(inSCE$scds_cxds_call) * 100)
        }

        if("scds_bcds_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "BCDS - Number of doublets",
                         "BCDS - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_bcds_call == TRUE),
                        sum(inSCE$scds_bcds_call == TRUE)/length(inSCE$scds_bcds_call) * 100)
        }

        if("scds_hybrid_call" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "SCDS Hybrid - Number of doublets",
                         "SCDS Hybrid - Percentage of doublets")
            values <- c(values, sum(inSCE$scds_hybrid_call == TRUE),
                        sum(inSCE$scds_hybrid_call == TRUE)/length(inSCE$scds_hybrid_call) * 100)
        }

        if("decontX_clusters" %in% colnames(SummarizedExperiment::colData(inSCE))){
            metrics <- c(metrics, "DecontX - Mean contamination",
                         "DecontX - Median contamination")
            values <- c(values, mean(inSCE$decontX_contamination),
                        stats::median(inSCE$decontX_contamination))
        }
    }

    df <- matrix(values)
    rownames(df) <- metrics
    colnames(df) <- colName
    return(df)
}

#' @title Plot table of SCTK QC outputs.
#' @description Plot QC metrics generated from QC algorithms via either kable or csv file.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' dimension reduction components or a variable with saved results. Required
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#'  'counts'.
#' @param simple Boolean. Indicates whether to generate a table of only
#' basic QC stats (ex. library size), or to generate a summary table of all
#' QC stats stored in the inSCE.
#' @param output Character vector. Options are "kable",
#' which a table will be generated via kable to an html file,
#' "csv", which a csv file will be generated, or "dataframe",
#' which outputs a dataframe.
#' @param fileName String. File name if output set to "csv".
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- sce[, colData(sce)$type != 'EmptyDroplet']
#' sampleSummaryStats(sce, simple = TRUE, output = "dataframe")
#' @importFrom magrittr %>%
#' @export
sampleSummaryStats <- function(inSCE,
                               sample = NULL,
                               useAssay = "counts",
                               simple = TRUE,
                               output = "kable",
                               fileName = "sctkSummaryStats.csv"){

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

    dfTableRes <- formatC(dfTableRes, drop0trailing = TRUE)

    if(output == "kable"){
        dfTableRes %>%
            knitr::kable(format = "html", align = rep('c', ncol(dfTableRes) + 1),
                         row.names = TRUE) %>% kableExtra::kable_styling()
    }else if(output == "csv"){
        utils::write.csv(dfTableRes, file = fileName)
    }else if(output == "dataframe"){
        return(dfTableRes)
    }

}
