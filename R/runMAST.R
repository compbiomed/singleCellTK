#' runMAST
#'
#' Run MAST for differential expression analysis on a
#' \linkS4class{SingleCellExperiment} inherited object.
#'
#' Condition group can be set by one and only one of the two ways: setting
#' \code{index1} and \code{index2}, or setting \code{class}, \code{classGroup1}
#' and \code{classGroup2}.
#'
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' which cells are being compared with those specified by \code{index1}. If
#' \code{NULL}, \code{index1} cells will be comapred with all other cells.
#' Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a single
#' string that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' will be compared with those specified by \code{classGroup1}. If \code{NULL},
#' \code{classGroup1} cells will be compared with all other cells.
#' @param comparisonName A string naming the DEG experiment. Will be required
#' in downstream plotting functions.
#' @param groupName1 A string naming the first comparison group - the condition
#' you are interested in.
#' @param groupName2 A string naming the second comparison group - the
#' condition that the first group is being compared with.
#' @param useThresh Whether to use adaptive thresholding to filter genes.
#' Default \code{FALSE}.
#' @param freqExpressed A numeric threshold that the genes expressed in less
#' than this fraction of cells will be removed. Default \code{0.1}.
#' @param onlyPos Whether to only output DEG with positive log2_FC value.
#' Default \code{FALSE}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this
#' value. Default \code{1}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$MAST} updated with the results: a list named by
#' \code{comparisonName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} storing the assay name that was used for
#' calculation, \code{$select} storing the cell selection indices (logical) for
#' each condition, and \code{$result} storing a \code{\link{data.frame}} of
#' the DEGs summary.
#' @export
#' @examples
#' data('sceBatches', package = 'singleCellTK')
#' sce.deg <- runMAST(sceBatches, class = "cell_type", comparisonName = 'aVSb',
#'                    classGroup1 = 'alpha', classGroup2 = 'beta',
#'                    groupName1 = 'a', groupName2 = 'b')
runMAST <- function(inSCE, useAssay = 'logcounts', index1 = NULL, index2 = NULL,
                    class = NULL, classGroup1 = NULL, classGroup2 = NULL,
                    comparisonName, groupName1 = NULL, groupName2 = NULL,
                    useThresh = FALSE, freqExpressed = 0.1, onlyPos = FALSE,
                    log2fcThreshold = NULL, fdrThreshold = 1){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop(paste('"useAssay" name: ', useAssay, ' not found.'))
    }
    if(is.null(index1) && (is.null(classGroup1) || is.null(class))){
        stop('At least "index1" or "classGroup1" should be specified.')
    } else if(!is.null(index1) && (!is.null(classGroup1) || !is.null(class))){
        stop('Only one of "index" and "class"/"classGroup1" ',
             'should be specified.')
    }
    if(!"MAST" %in% names(S4Vectors::metadata(inSCE))){
        S4Vectors::metadata(inSCE)$MAST <- list()
    }
    if(comparisonName %in% names(S4Vectors::metadata(inSCE)$MAST)){
        warning(paste0('comparisonName: "', comparisonName, '", already '),
                'exists. Overwriting.')
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
    if(all(is.na(SummarizedExperiment::colData(sca)$cngeneson))){
        SummarizedExperiment::colData(sca)$cngeneson <- 0
    }
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
    fcHurdleSig <- merge(fcHurdle,
                         data.table::as.data.table(S4Vectors::mcols(sca)),
                         by = "primerid")
    if(!is.null(log2fcThreshold)){
        fcHurdleSig <- fcHurdleSig[abs(fcHurdleSig$coef) > log2fcThreshold,]
    }
    if(isTRUE(onlyPos)){
        fcHurdleSig <- fcHurdleSig[fcHurdleSig$coef > 0,]
    }
    if(!is.null(fdrThreshold)){
        fcHurdleSig <- fcHurdleSig[fcHurdleSig$fdr < fdrThreshold,]
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
             result = data.frame(fcHurdleSig))
    return(inSCE)
}

#' Find the marker gene set for each cluster
#' With an input SingleCellExperiment object and specifying the clustering
#' labels, this function iteratively call the differential expression analysis
#' on each cluster against all the others.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' MAST calculations. Default \code{"logcounts"}.
#' @param cluster One single character to specify a column in
#' \code{colData(inSCE)} for the clustering label. Alternatively, a vector or
#' a factor is also acceptable. Default \code{"cluster"}.
#' @param useThresh Whether to use adaptive thresholding to filter genes.
#' Default \code{FALSE}.
#' @param freqExpressed A numeric threshold that the genes expressed in less
#' than this fraction of cells will be removed. Default \code{0.1}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}
#' @param fdrThreshold Only out put DEGs with FDR value smaller than this
#' value. Default \code{1}
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$findMarker} updated with a data.table of the up-
#' regulated DEGs for each cluster.
#' @author Yichen Wang
findMarkerDiffExp <- function(inSCE, useAssay = 'logcounts',
                              cluster = 'cluster', log2fcThreshold = NULL,
                              fdrThreshold = 1, useThresh = FALSE,
                              freqExpressed = 0.1){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop('"useAssay" name: ', useAssay, ' not found.')
    }
    ## Check whether use colData or customized vector/factor
    if(is.character(cluster) && length(cluster) == 1){
        if(!cluster %in% names(SummarizedExperiment::colData(inSCE))){
            stop('"cluster": ', cluster, ' not found in colData(inSCE).')
        }
        clusterName <- cluster
        cluster <- SummarizedExperiment::colData(inSCE)[[cluster]]
    } else if(!length(cluster) == ncol(inSCE)){
        stop('The number of cluster labels does not match the number of cells.')
    } else {
        clusterName <- 'findMarker_cluster'
        SummarizedExperiment::colData(inSCE)[[clusterName]] <- cluster
    }
    # Iterate
    if(is.factor(cluster)){
        # In case inSCE is a subset, when "levels" is a full list of all
        # clusters, want to keep the ordering stored in the factor.
        uniqClust <- as.vector(unique(cluster))
        levelOrder <- levels(cluster)
        uniqClust <- levelOrder[levelOrder %in% uniqClust]
    } else {
        uniqClust <- unique(cluster)
    }
    for(c in uniqClust){
        message('Computing for cluster: ', c)
        clusterIndex <- cluster == c
        inSCE <- runMAST(inSCE, useAssay = useAssay, index1 = clusterIndex,
                         comparisonName = paste0('findMarker', c),
                         onlyPos = TRUE, log2fcThreshold = log2fcThreshold,
                         fdrThreshold = fdrThreshold, useThresh = useThresh,
                         freqExpressed = freqExpressed)
    }
    degFull <- NULL
    for(c in uniqClust){
        degTable <-
            S4Vectors::metadata(inSCE)$MAST[[paste0('findMarker', c)]]$result
        degTable[[clusterName]] <- c
        if(is.null(degFull)){
            degFull <- degTable
        } else {
            degFull <- rbind(degFull, degTable)
        }
    }
    S4Vectors::metadata(inSCE)$MAST[paste0('findMarker', uniqClust)] <- NULL
    if(length(names(S4Vectors::metadata(inSCE)$MAST)) == 0){
        S4Vectors::metadata(inSCE)$MAST <- NULL
    }
    S4Vectors::metadata(inSCE)$findMarker <- degFull
    return(inSCE)
}
