#' label cell types with SingleR
#' @description to be filled in
#' @param inSCE In put SCE object
#' @param useAssay The assay to calculate
#' @param useRef Specify an reference provided by SingleR
#' @param level Cell type labeling level
#' @param featureType Gene symbols or Ensembl IDs
#' @return An SCE object
#' @export
runSingleR <- function(inSCE, useAssay = "logcounts",
                       useRef = c("hpca", "be", "mp", "dice",
                                  "nh", "mi", "ig", "mrs"),
                       level = c("main", "fine", "ont"),
                       featureType = c("symbol", "ensembl")) {
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop('"useAssay" name: ', useAssay, ' not found.')
    }

    useRef <- match.arg(useRef)
    featureType <- match.arg(featureType)
    level <- match.arg(level)

    if (featureType == "symbol") {
        useEnsembl <- FALSE
    } else {
        useEnsembl <- TRUE
    }

    if (useRef == "hpca") {
        message("Loading reference data 'HumanPrimaryCellAtlasData'...")
        ref <- SingleR::HumanPrimaryCellAtlasData(ensembl = useEnsembl,
                                                  cell.ont = "none")
        labelColName <- paste0("label.", level)
    } else if (useRef == "be") {
        message("Loading reference data 'BlueprintEncodeData'...")
        ref <- SingleR::BlueprintEncodeData(ensembl = useEnsembl,
                                            cell.ont = "none")
        labelColName <- paste0("label.", level)
    } else if (useRef == "mp") {
        message("Loading reference data 'MuraroPancreasData'...")
        ref <- scRNAseq::MuraroPancreasData(ensembl = useEnsembl)
        if (!isTRUE(useEnsembl)) {
            rownames(ref) <- SummarizedExperiment::rowData(ref)$symbol
        }
        ref <- ref[,!is.na(ref$label) & ref$label!="unclear"]
        ref <- scater_logNormCounts(ref, logAssayName = "logcounts")
        labelColName <- "label"
        warning("MuraroPancreasData does not have multiple levels of label. ",
                "Using its default labeling.")
    }

    predictions <- SingleR::SingleR(test = inSCE, assay.type.test = useAssay,
                                    ref = ref,
                                    labels = ref[[labelColName]])
    colnames(predictions) <- paste0("SingleR_", useRef, "_", level, "_",
                                    colnames(predictions))
    SummarizedExperiment::colData(inSCE) <-
        cbind(SummarizedExperiment::colData(inSCE), predictions)

    return(inSCE)
}
