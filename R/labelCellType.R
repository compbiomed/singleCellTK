#' label cell types with SingleR
#' @description to be filled in
#' @param inSCE In put SCE object
#' @param useAssay The assay to calculate
#' @param useSCERef A customized reference dataset in SCE format
#' @param labelColName The colData column for cell type labeling, stored in
#' useSCERef
#' @param useBltinRef Specify an reference provided by SingleR
#' @param level Cell type labeling level, only supported by a few SingleR
#' references.
#' @param featureType Use gene symbols or Ensembl IDs, only supported by SingleR
#' references.
#' @param labelByCluster Label inSCE clusters already identified, instead of
#' labeling each cells.
#' @return An SCE object
#' @export
runSingleR <- function(inSCE,
                       useAssay = "logcounts",
                       useSCERef = NULL,
                       labelColName = NULL,
                       useBltinRef = c("hpca", "bpe", "mp", "dice",
                                  "immgen", "mouse", "zeisel"),
                       level = c("main", "fine", "ont"),
                       featureType = c("symbol", "ensembl"),
                       labelByCluster = NULL) {
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop('"useAssay" name: ', useAssay, ' not found.')
    }

    useBltinRef <- match.arg(useBltinRef)
    featureType <- match.arg(featureType)
    level <- match.arg(level)

    if (featureType == "symbol") {
        useEnsembl <- FALSE
    } else {
        useEnsembl <- TRUE
    }

    if (!is.null(useSCERef)) {
        ref <- useSCERef
        if (is.null(labelColName)) {
            stop("labelColName must be specified if given an SCE reference.")
        } else {
            if (!labelColName %in% SummarizedExperiment::colData(ref)) {
                stop("Specified labelColName not found in given reference.")
            }
            message("Customized reference does not support level setting. ",
                    "Please check if the colData using matches your needs.")
            message("Ensembl vs symbol transfering not supported for cusomized",
                    " reference. Please make sure its rownames match with ",
                    "inSCE.")
        }
    } else {
        if (useBltinRef == "hpca") {
            message("Loading reference data 'HumanPrimaryCellAtlasData'...")
            ref <- SingleR::HumanPrimaryCellAtlasData(ensembl = useEnsembl,
                                                      cell.ont = "none")
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "bpe") {
            message("Loading reference data 'BlueprintEncodeData'...")
            ref <- SingleR::BlueprintEncodeData(ensembl = useEnsembl,
                                                cell.ont = "none")
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "mp") {
            message("Loading reference data 'MuraroPancreasData'...")
            ref <- scRNAseq::MuraroPancreasData(ensembl = useEnsembl)
            if (!isTRUE(useEnsembl)) {
                rownames(ref) <- SummarizedExperiment::rowData(ref)$symbol
            }
            ref <- ref[,!is.na(ref$label) & ref$label!="unclear"]
            ref <- scater_logNormCounts(ref, logAssayName = "logcounts")
            labelColName <- "label"
            warning("MuraroPancreasData does not have multiple levels of ",
                    "label. Using its default labeling.")
        } else if (useBltinRef == "dice") {
            ref <- celldex::DatabaseImmuneCellExpressionData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "immgen") {
            ref <- celldex::ImmGenData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "mouse") {
            ref <- celldex::MouseRNAseqData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "zeisel") {
            message("Loading reference data 'MuraroPancreasData'...")
            ref <- scRNAseq::ZeiselBrainData(ensembl = useEnsembl)
            ref <- ref[,ref$level2class!="(none)"]
            ref <- scater_logNormCounts(ref, logAssayName = "logcounts")
            labelColName <- "level2class"
            warning("ZeiselBrainData does not support levels. ",
                    "Using its default labeling.")
        }
    }
    if (is.null(labelByCluster)) {
        clusters <- NULL
    } else {
        clusters <- inSCE[[labelByCluster]]
    }
    predictions <- SingleR::SingleR(test = inSCE, assay.type.test = useAssay,
                                    ref = ref, clusters = clusters,
                                    labels = ref[[labelColName]])
    predictions$tuning.scores <- NULL
    if (is.null(clusters)) {
        colnames(predictions) <- paste0("SingleR_", useBltinRef, "_", level,
                                        "_", colnames(predictions))
        for (n in colnames(predictions)) {
            inSCE[[n]] <- predictions[[n]]
        }
    } else {
        colnames(predictions) <- paste0("SingleR_", useBltinRef, "_", level,
                                        "_", colnames(predictions))
        for (i in seq(nrow(predictions))) {
            clusterLabel <- rownames(predictions)[i]
            for (n in colnames(predictions)) {
                inSCE[[n]][clusters == clusterLabel] <- predictions[[n]][i]
            }
        }
    }

    return(inSCE)
}

