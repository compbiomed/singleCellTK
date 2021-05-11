#' Label cell types with SingleR
#' @description
#' SingleR works with a reference dataset where the cell type
#' labeling is given. Given a reference dataset of samples (single-cell or bulk)
#' with known labels, it assigns those labels to new cells from a test dataset
#' based on similarities in their expression profiles.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay character. A string specifying which assay to use for
#' expression profile identification. Required.
#' @param useSCERef \linkS4class{SingleCellExperiment} inherited object. An
#' optional customized reference dataset. Default \code{NULL}.
#' @param labelColName A single character. A string specifying the column in
#' \code{colData(useSCERef)} that stores the cell type labeling. Default
#' \code{NULL}.
#' @param useBltinRef A single character. A string that specifies a reference
#' provided by SingleR. Choose from \code{"hpca", "bpe", "mp", "dice", "immgen",
#' "mouse", "zeisel"}. See detail. Default \code{"hpca"}.
#' @param level A string for cell type labeling level. Used only when using
#' some of the SingleR built-in references. Choose from \code{"main", "fine",
#' "ont"}. Default \code{"main"}.
#' @param featureType A string for whether to use gene symbols or Ensembl IDs
#' when using a SingleR built-in reference. Should be set based on the type of
#' \code{rownames} of \code{inSCE}. Choose from \code{"symbol", "ensembl"}.
#' Default \code{"symbol"}.
#' @param labelByCluster A single character. A string specifying the column name
#' in \code{colData(inSCE)} that stores clustering labels. Use this when users
#' want to only label cells on cluster level, instead of performing calculation
#' on each cell. Default \code{NULL}.
#' @return Input SCE object with cell type labeling updated in
#' \code{colData(inSCE)}, together with scoring metrics.
#' @export
#' @examples
#' data("sceBatches")
#' logcounts(sceBatches) <- log(counts(sceBatches) + 1)
#' #sceBatches <- runSingleR(sceBatches, useBltinRef = "mp")
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
    if(!useAssay %in% expDataNames(inSCE)){
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
            ref <- celldex::HumanPrimaryCellAtlasData(ensembl = useEnsembl,
                                                      cell.ont = "none")
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "bpe") {
            message("Loading reference data 'BlueprintEncodeData'...")
            ref <- celldex::BlueprintEncodeData(ensembl = useEnsembl,
                                                cell.ont = "none")
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "mp") {
            message("Loading reference data 'MuraroPancreasData'...")
            ref <- scRNAseq::MuraroPancreasData(ensembl = useEnsembl)
            if (!isTRUE(useEnsembl)) {
                rownames(ref) <- SummarizedExperiment::rowData(ref)$symbol
            }
            ref <- ref[,!is.na(ref$label) & ref$label!="unclear"]
            ref <- scaterlogNormCounts(ref, assayName = "logcounts")
            labelColName <- "label"
            warning("MuraroPancreasData does not have multiple levels of ",
                    "label. Using its default labeling.")
        } else if (useBltinRef == "dice") {
            message("Loading reference data 'DatabaseImmuneCellExpressionData'...")
            ref <- celldex::DatabaseImmuneCellExpressionData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "immgen") {
            message("Loading reference data 'ImmGenData'...")
            ref <- celldex::ImmGenData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "mouse") {
            message("Loading reference data 'MouseRNAseqData'...")
            ref <- celldex::MouseRNAseqData(ensembl = useEnsembl)
            labelColName <- paste0("label.", level)
        } else if (useBltinRef == "zeisel") {
            message("Loading reference data 'ZeiselBrainData'...")
            ref <- scRNAseq::ZeiselBrainData(ensembl = useEnsembl)
            ref <- ref[,ref$level2class!="(none)"]
            ref <- scaterlogNormCounts(ref, assayName = "logcounts")
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
    # predictions <- SingleR::SingleR(test = inSCE, assay.type.test = useAssay,
    #                                 ref = ref, clusters = clusters,
    #                                 labels = ref[[labelColName]])
    predictions <- SingleR::SingleR(test = expData(inSCE, useAssay),
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

