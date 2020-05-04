#' runMAST
#'
#' Run MAST for differential expression analysis on a SCtkExperiment object.
#'
#' @param inSCE SCtkExperiment object
#' @param useAssay character, default "logcounts". A string specifying which
#' assay to use for the MAST calculations.
#' @param index1 character/numeric/logical, default NULL. A vector that
#' specifies which cells are of interests.
#' @param index2 character/numeric/logical, default NULL. A vector that
#' specifies which cells are being compared with those specified by `index1`.
#' If NULL, `index1` cells will be comapred with all other cells.
#' @param class character/numeric or factor, default NULL. A vector/factor of
#' `ncol(inSCE)` elements, or a single string that specifies a colname of
#' `colData(inSCE)`.
#' @param classGroup1 character/numeric, default NULL. Specifying which
#' "levels" given in `class` are of interests.
#' @param classGroup2 character/numeric, default NULL. Specifying which
#' "levels" given in `class` will be compared with those specified by
#' `classGroup1`. If NULL, `classGroup1` cells will be compared with all other
#' cells.
#' @param comparisonName character. A string naming the DEG experiment.
#' @param groupName1 character. A string naming the first comparison group.
#' @param groupName2 character. A string naming the second comparison group.
#' @param useThresh logical, default FALSE. Whether to use adaptive
#' thresholding to filter genes.
#' @param freqExpressed numeric, default 0.1. A threshold that the genes
#' expressed in less than this fraction of cells will be removed.
#' @param onlyPos logical, default FALSE. Whether to only output DEG with
#' positive log2_FC value.
#' @param log2fcThreshold numeric, default NULL. Only out put DEGs with the
#' absolute values of log2FC larger than this value.
#' @param fdrThreshold numeric, default 1. Only out put DEGs with FDR value
#' smaller than this value.
#' @return The input SCtkExperiment object with `metadata(inSCE)$MAST`` updated
#' with the results: a list named by `comparisonName`, with `groupNames`
#' containing the naming of the two conditions, `useAssay` storing the assay
#' name that was used for calculation, `select` storing the cell selection
#' indices (logical) for each condition, and `result` storing a `data.frame` of
#' the DEGs summary.
#' @export
#' @examples
#' data(sceBatch)
#' sce.deg <- runMAST(sceBatch, class = "cell_type", comparisonName = 'aVSb',
#'                    classGroup1 = 'alpha', classGroup2 = 'beta')
runMAST <- function(inSCE, useAssay = 'logcounts', index1 = NULL, index2 = NULL,
                    class = NULL, classGroup1 = NULL, classGroup2 = NULL,
                    comparisonName, groupName1, groupName2, useThresh = FALSE,
                    freqExpressed = 0.1, onlyPos = FALSE,
                    log2fcThreshold = NULL, fdrThreshold = 1){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop("\"inSCE\" should be a SCtkExperiment Object.")
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop(paste("\"useAssay\" name: ", useAssay, " not found."))
    }
    if(is.null(index1) && (is.null(classGroup1) || is.null(class))){
        stop("At least \"index1\" or \"classGroup1\" should be specified.")
    } else if(!is.null(index1) && (!is.null(classGroup1) || !is.null(class))){
        stop("Only one of \"index\" and \"class\"/\"classGroup1\" ",
             "should be specified.")
    }
    if(!"MAST" %in% names(S4Vectors::metadata(inSCE))){
        S4Vectors::metadata(inSCE)$MAST <- list()
    }
    if(comparisonName %in% names(S4Vectors::metadata(inSCE)$MAST)){
        warning(paste0("comparisonName: \"", comparisonName, "\", already "),
                "exists. Overwriting.")
    }
    groupNames <- c(groupName1, groupName2)

    # Data preparation
    ## Identify the subsetting
    if(!is.null(index1)){
        cells1 <- colnames(inSCE[,index1])
        if(!is.null(index2)){
            cells2 <- colnames(inSCE[,index2])
        } else {
            cells2 <- sort(setdiff(colnames(inSCE), cells1))
        }
    } else {
        if(length(class) == 1 && class(class) == "character"){
            class <- SummarizedExperiment::colData(inSCE)[[class]]
        }
        index1 <- class %in% classGroup1
        if(is.null(classGroup2)){
            index2 <- !class %in% classGroup1
        } else {
            index2 <- class %in% classGroup2
        }
        cells1 <- colnames(inSCE[,index1])
        cells2 <- colnames(inSCE[,index2])
    }
    if(length(cells1) == 0){
      stop("Number of cells selected for group1 equals to zero.")
    }
    if(length(cells2) == 0){
      stop("Number of cells selected for group2 equals to zero.")
    }
    ix1 <- colnames(inSCE) %in% cells1
    ix2 <- colnames(inSCE) %in% cells2
    select <- list(ix1 = ix1, ix2 = ix2)

    ## Extract
    mat <-
        SummarizedExperiment::assay(inSCE, useAssay)[,c(cells1, cells2)]
    if(!inherits(mat, 'matrix')){
        mat <- as.matrix(mat)
    }
    mat <- featureNameDedup(mat)
    cdat <- data.frame(wellKey = colnames(mat),
        condition = as.factor(c(rep("c1", length(cells1)),
                                rep("c2", length(cells2)))),
        ngeneson = rep("", (length(cells1) + length(cells2))),
        stringsAsFactors = FALSE)
    sca <- MAST::FromMatrix(mat, cdat)
    cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
    SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
    cond <- factor(SummarizedExperiment::colData(sca)$condition)
    cond <- stats::relevel(cond, "c2")
    SummarizedExperiment::colData(sca)$condition <- cond
    # Calculation
    ## Filtration
    if(useThresh){
        sca <- sca[which(MAST::freq(sca) > 0),]
        invisible(utils::capture.output(thresh <-
            MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(sca),
                                            nbins = 20, min_per_bin = 30)))
        SummarizedExperiment::assays(sca) <-
            list(thresh = thresh$counts_threshold,
                 tpm = SummarizedExperiment::assay(sca))
    }
    if(sum(MAST::freq(sca) > freqExpressed) <= 1){
        stop("Not enough genes pass frequency expressed filter of 1")
    }
    sca <- sca[which(MAST::freq(sca) > freqExpressed),]
    message(paste("Using", nrow(sca), 'genes after filtration.'))
    ##
    SummarizedExperiment::colData(sca)$cngeneson <-
        scale(colSums(SummarizedExperiment::assay(sca) > 0))
    zlmCond <- MAST::zlm(~condition + cngeneson, sca)
    summaryCond <- MAST::summary(zlmCond, doLRT = "conditionc1")
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[summaryDt$contrast == "conditionc1" &
                                summaryDt$component == "H", 
                                c('primerid', 'Pr(>Chisq)')],
                      summaryDt[summaryDt$contrast == "conditionc1" & 
                                summaryDt$component == "logFC", 
                                c('primerid', 'coef', 'ci.hi', 'ci.lo')],
                      by = "primerid")
    fcHurdle$fdr <- stats::p.adjust(fcHurdle$`Pr(>Chisq)`, "fdr")
    if (is.null(log2fcThreshold)) {
      fcHurdleSig <- fcHurdle
    } else {
      fcHurdleSig <- merge(fcHurdle[fcHurdle$fdr < fdrThreshold & 
                                    abs(fcHurdle$coef) > log2fcThreshold], 
                                    data.table::as.data.table(S4Vectors::mcols(sca)),
                           by = "primerid")
      if (onlyPos) {
        fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc >
            0), ]
      }
    }
    fcHurdleSig <- fcHurdleSig[, -c(4, 5)]
    names(fcHurdleSig)[c(1, 2, 3, 4)] <- c("Gene", "Pvalue",
      "Log2_FC", "FDR")
    fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$Pvalue, decreasing = FALSE),
      ]
    S4Vectors::metadata(inSCE)$MAST[[comparisonName]] <-
        list(groupNames = groupNames,
             useAssay = useAssay,
             select = select,
             result = fcHurdleSig)
    return(inSCE)
}
