#' Helper function for differential expression analysis methods that accepts
#' multiple ways of conditional subsetting and returns stable index format.
#' Meanwhile it does all the input checkings.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay character. A string specifying which assay to use. Required.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' the control group against those specified by \code{index1}. If
#' \code{NULL} when using index specification, \code{index1} cells will be
#' compared with all other cells. Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a character
#' scalar that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' is the control group against those specified by \code{classGroup1}. If
#' \code{NULL} when using annotation specification, \code{classGroup1} cells
#' will be compared with all other cells.
#' @param groupName1 A character scalar naming the group of interests. Required.
#' @param groupName2 A character scalar naming the control group. Required.
#' @param analysisName A character scalar naming the DEG analysis. Required
#' @param covariates A character vector of additional covariates used in linear
#' regression methods such as Limma and DESeq2. Default \code{NULL}
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @return A list object with part of formatted DE analysis information
#' @author Yichen Wang
.formatDEAList <- function(inSCE, useAssay, index1 = NULL, index2 = NULL,
                           class = NULL, classGroup1 = NULL,
                           classGroup2 = NULL, groupName1, groupName2,
                           analysisName, covariates = NULL, overwrite = FALSE){
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
    if(!is.null(covariates) &&
       !all(covariates %in% names(SummarizedExperiment::colData(inSCE)))){
        stop("Not all specified covariates exist.")
    }
    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
        if(analysisName %in% names(S4Vectors::metadata(inSCE)$diffExp)){
            if(isTRUE(overwrite)){
                warning("analysisName '", analysisName, "' already exists, ",
                        "will be overwritten.")
            } else {
                stop("analysisName '", analysisName, "' already exists. ",
                     "Set `overwrite` to `TRUE` to overwrite.")
            }
        }
    }
    groupNames <- c(groupName1, groupName2)
    annotation <- c(class, covariates)
    if(!is.null(index1)){
        cells1 <- colnames(inSCE[,index1])
        if(!is.null(index2)){
            cells2 <- colnames(inSCE[,index2])
        } else {
            cells2 <- sort(setdiff(colnames(inSCE), cells1))
        }
    } else {
        if(length(class) == 1 && inherits(class, "character")){
            if(!class %in% names(SummarizedExperiment::colData(inSCE))){
                stop("class: '", class, "' not found.")
            }
            class <- SummarizedExperiment::colData(inSCE)[[class]]
        } else {
            if(!length(class) == ncol(inSCE)){
                stop("Length of given `class` vector should equal to ",
                     "`ncol(inSCE)`; Or specify a column name of ",
                     "`colData(inSCE)`")
            }
        }
        uniqCats <- unique(as.vector(class))
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
    if(!is.null(covariates)){
        for (c in covariates){
            col <- SummarizedExperiment::colData(inSCE)[(ix1 | ix2),c]
            if(length(unique(as.vector(col))) < 2){
                stop("Less than 2 levels in specified covariate: ", c)
            }
        }
    }
    return(
        list(useAssay = useAssay,
             groupNames = groupNames,
             select = select,
             annotation = annotation)
    )
}

#' Perform differential expression analysis on SCE with specified method
#' Method supported: 'MAST', 'DESeq2', 'Limma', 'ANOVA'
#' @param method A single character for specific method. Choose from 'MAST',
#' 'DESeq2', 'Limma', 'ANOVA'. Required
#' @param ... Other arguments passed to specific functions. Refer to
#' \code{\link{runMAST}}, \code{\link{runDESeq2}}, \code{\link{runLimmaDE}},
#' \code{\link{runANOVA}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDEAnalysis(inSCE = sce, groupName1 = "Sample1", method = "DESeq2",
#'  groupName2 = "Sample2", index1 = 1:20, index2 = 21:40,
#'  analysisName = "DESeq2")
#' @return Input SCE object with \code{metadata(inSCE)} updated with name
#' \code{"diffExp"} as a \code{list} object. Detail refers to the four child
#' functions.
#' @export
runDEAnalysis <- function(method = 'MAST', ...){
    if(!method %in% c('MAST', 'DESeq2', 'Limma', 'ANOVA')){
        stop("method should be one of: 'MAST', 'DESeq2', 'Limma', 'ANOVA'")
    }
    funcList <- list(MAST = runMAST,
                     DESeq2 = runDESeq2,
                     Limma = runLimmaDE,
                     ANOVA = runANOVA)
    funcList[[method]](...)
}

#' Perform differential expression analysis on SCE with DESeq2.
#'
#' Condition specification allows two methods:
#' 1. Index level selection. Arguments \code{index1} and \code{index2} will be
#' used.
#' 2. Annotation level selection. Arguments \code{class}, \code{classGroup1} and
#' \code{classGroup2} will be used.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' DESeq2 regression. Default \code{"logcounts"}.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' the control group against those specified by \code{index1}. If
#' \code{NULL} when using index specification, \code{index1} cells will be
#' compared with all other cells. Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a character
#' scalar that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' is the control group against those specified by \code{classGroup1}. If
#' \code{NULL} when using annotation specification, \code{classGroup1} cells
#' will be compared with all other cells.
#' @param analysisName A character scalar naming the DEG analysis. Required
#' @param groupName1 A character scalar naming the group of interests. Required.
#' @param groupName2 A character scalar naming the control group. Required.
#' @param covariates A character vector of additional covariates to use when
#' building the model. All covariates must exist in
#' \code{names(colData(inSCE))}. Default \code{NULL}.
#' @param fullReduced Whether to apply LRT (Likelihood ratio test) with a 'full'
#' model. Default \code{TRUE}.
#' @param onlyPos Whether to only output DEG with positive log2_FC value.
#' Default \code{FALSE}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' greater than this value. Default \code{0.25}
#' @param fdrThreshold Only out put DEGs with FDR value less than this
#' value. Default \code{0.05}
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDESeq2(inSCE = sce, groupName1 = "Sample1",
#'  groupName2 = "Sample2", index1 = 1:20, index2 = 21:40,
#'  analysisName = "DESeq2")
#'
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$DESeq2} updated with the results: a list named by
#' \code{analysisName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} storing the assay name that was used for
#' calculation, \code{$select} storing the cell selection indices (logical) for
#' each condition, \code{$result} storing a \code{\link{data.frame}} of
#' the DEGs summary, and \code{$method} storing \code{"DESeq2"}.
#' @export
runDESeq2 <- function(inSCE, useAssay = 'counts', index1 = NULL,
                      index2 = NULL, class = NULL, classGroup1 = NULL,
                      classGroup2 = NULL, analysisName, groupName1,
                      groupName2, covariates = NULL, fullReduced = TRUE,
                      onlyPos = FALSE, log2fcThreshold = NULL,
                      fdrThreshold = 1, overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, index1, index2, class,
                                 classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)
    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    subsetIdx <- (ix1 | ix2)
    conditions <- rep(NA, ncol(inSCE))
    conditions[ix1] <- 'cond1'
    conditions[ix2] <- 'cond2'
    conditions <- conditions[!is.na(conditions)]
    annotData <- data.frame(condition = conditions,
                            row.names = colnames(inSCE)[subsetIdx])
    cov <- SummarizedExperiment::colData(inSCE)[subsetIdx, covariates,
                                                drop = FALSE]
    annotData <- cbind(annotData, cov)
    mat <- SummarizedExperiment::assay(inSCE[,subsetIdx], useAssay)
    if(!inherits(mat, 'matrix')){
        mat <- as.matrix(mat)
    }
    mat <- featureNameDedup(mat)
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = mat, colData = annotData,
        design = stats::as.formula(paste0("~", c('condition', covariates),
                                          collapse = "+"))
    )
    if(isTRUE(fullReduced)){
        dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1)
    } else {
        dds <- DESeq2::DESeq(
            dds, test = "LRT",
            reduced = stats::as.formula(paste0('~condition', covariates,
                                               collapse = '+'))
        )
    }
    res <- DESeq2::results(dds, pAdjustMethod = 'fdr')
    deg <- data.frame(res)[,c(-1, -3, -4)]
    deg <- cbind(data.frame(Gene = rownames(deg), stringsAsFactors = FALSE),
                 deg)
    rownames(deg) <- NULL
    colnames(deg) <- c('Gene', 'Log2_FC', 'Pvalue', 'FDR')
    deg$Log2_FC <- - deg$Log2_FC # JNMLP
    # Result Filtration
    if(isTRUE(onlyPos)){
        deg <- deg[deg$Log2_FC > 0,]
    }
    if(!is.null(fdrThreshold)){
        deg <- deg[deg$FDR < fdrThreshold,]
    }
    if(!is.null(log2fcThreshold)){
        deg <- deg[abs(deg$Log2_FC) > log2fcThreshold,]
    }
    # Format output
    deg <- deg[order(deg$FDR, na.last = TRUE),]
    resultList$result <- deg
    resultList$method <- 'DESeq2'
    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    } else {
        S4Vectors::metadata(inSCE)$diffExp <- list()
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    }
    return(inSCE)
}


#' Perform differential expression analysis on SCE with Limma.
#'
#' Condition specification allows two methods:
#' 1. Index level selection. Arguments \code{index1} and \code{index2} will be
#' used.
#' 2. Annotation level selection. Arguments \code{class}, \code{classGroup1} and
#' \code{classGroup2} will be used.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for the
#' Limma regression. Default \code{"logcounts"}.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' the control group against those specified by \code{index1}. If
#' \code{NULL} when using index specification, \code{index1} cells will be
#' compared with all other cells. Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a character
#' scalar that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' is the control group against those specified by \code{classGroup1}. If
#' \code{NULL} when using annotation specification, \code{classGroup1} cells
#' will be compared with all other cells.
#' @param analysisName A character scalar naming the DEG analysis. Required
#' @param groupName1 A character scalar naming the group of interests. Required.
#' @param groupName2 A character scalar naming the control group. Required.
#' @param covariates A character vector of additional covariates to use when
#' building the model. All covariates must exist in
#' \code{names(colData(inSCE))}. Default \code{NULL}.
#' @param onlyPos Whether to only output DEG with positive log2_FC value.
#' Default \code{FALSE}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' greater than this value. Default \code{0.25}
#' @param fdrThreshold Only out put DEGs with FDR value less than this
#' value. Default \code{0.05}
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' library(scater)
#' sce <- logNormCounts(sce)
#' sce <- runLimmaDE(inSCE = sce, groupName1 = "Sample1",
#'  groupName2 = "Sample2", index1 = 1:20, index2 = 21:40,
#'  analysisName = "Limma")
#'
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$diffExp} updated with the results: a list named by
#' \code{analysisName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} storing the assay name that was used for
#' calculation, \code{$select} storing the cell selection indices (logical) for
#' each condition, \code{$result} storing a \code{\link{data.frame}} of
#' the DEGs summary, and \code{$method} storing \code{"Limma"}.
#' @export
runLimmaDE <- function(inSCE, useAssay = 'logcounts', index1 = NULL,
                       index2 = NULL, class = NULL, classGroup1 = NULL,
                       classGroup2 = NULL, analysisName, groupName1,
                       groupName2, covariates = NULL, onlyPos = FALSE,
                       log2fcThreshold = 0.25, fdrThreshold = 0.05,
                       overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, index1, index2, class,
                                 classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)
    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    subsetIdx <- (ix1 | ix2)
    conditions <- rep(NA, ncol(inSCE))
    conditions[ix1] <- 'cond1'
    conditions[ix2] <- 'cond2'
    conditions <- conditions[!is.na(conditions)]
    annotData <- data.frame(condition = conditions,
                            row.names = colnames(inSCE)[subsetIdx])
    cov <- SummarizedExperiment::colData(inSCE)[subsetIdx, covariates,
                                                drop = FALSE]
    annotData <- cbind(annotData, cov)
    mat <- SummarizedExperiment::assay(inSCE[,subsetIdx], useAssay)
    if(!inherits(mat, 'matrix')){
        mat <- as.matrix(mat)
    }
    mat <- featureNameDedup(mat)
    design <- stats::model.matrix(
        stats::as.formula(paste0("~", paste0(c('condition', covariates),
                                             collapse = "+"))),
        data = annotData)
    fit <- limma::lmFit(mat, design)
    ebayes <- limma::eBayes(fit)
    deg <- limma::topTable(ebayes, adjust = 'fdr',
                           number = nrow(inSCE))
    deg <- deg[,c(-2, -3, -6)]
    deg <- cbind(data.frame(Gene = rownames(deg), stringsAsFactors = FALSE),
                 deg)
    rownames(deg) <- NULL
    colnames(deg) <- c("Gene", "Log2_FC", "Pvalue", "FDR")
    deg$Log2_FC <- - deg$Log2_FC # JNMLP
    # Result Filtration
    if(isTRUE(onlyPos)){
      deg <- deg[deg$Log2_FC > 0,]
    }
    if(!is.null(fdrThreshold)){
      deg <- deg[deg$FDR < fdrThreshold,]
    }
    if(!is.null(log2fcThreshold)){
      deg <- deg[abs(deg$Log2_FC) > log2fcThreshold,]
    }
    # Format output
    deg <- deg[order(deg$FDR, na.last = TRUE),]
    resultList$result <- deg
    resultList$method <- 'Limma'
    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
      S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    } else {
      S4Vectors::metadata(inSCE)$diffExp <- list()
      S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    }
    return(inSCE)
}

#' Perform differential expression analysis on SCE with ANOVA
#'
#' Condition specification allows two methods:
#' 1. Index level selection. Arguments \code{index1} and \code{index2} will be
#' used.
#' 2. Annotation level selection. Arguments \code{class}, \code{classGroup1} and
#' \code{classGroup2} will be used.
#'
#' NOTE that ANOVA method does not produce Log2FC value, but P-value and FDR
#' only.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for ANOVA.
#' Default \code{"logcounts"}.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' the control group against those specified by \code{index1}. If
#' \code{NULL} when using index specification, \code{index1} cells will be
#' compared with all other cells. Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a character
#' scalar that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' is the control group against those specified by \code{classGroup1}. If
#' \code{NULL} when using annotation specification, \code{classGroup1} cells
#' will be compared with all other cells.
#' @param analysisName A character scalar naming the DEG analysis. Required
#' @param groupName1 A character scalar naming the group of interests. Required.
#' @param groupName2 A character scalar naming the control group. Required.
#' @param covariates A character vector of additional covariates to use when
#' building the model. All covariates must exist in
#' \code{names(colData(inSCE))}. Default \code{NULL}.
#' @param onlyPos Whether to only output DEG with positive log2_FC value.
#' Default \code{FALSE}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' greater than this value. Default \code{0.25}
#' @param fdrThreshold Only out put DEGs with FDR value less than this
#' value. Default \code{0.05}
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' library(scater)
#' sce <- logNormCounts(sce)
#' sce <- runANOVA(inSCE = sce, groupName1 = "Sample1",
#'  groupName2 = "Sample2", index1 = 1:20, index2 = 21:40,
#'  analysisName = "ANOVA", fdrThreshold = NULL)
#'
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$diffExp} updated with the results: a list named by
#' \code{analysisName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} storing the assay name that was used for
#' calculation, \code{$select} storing the cell selection indices (logical) for
#' each condition, \code{$result} storing a \code{\link{data.frame}} of
#' the DEGs summary, and \code{$method} storing \code{"ANOVA"}.
#' @export
runANOVA <- function(inSCE, useAssay = 'logcounts', index1 = NULL,
                     index2 = NULL, class = NULL, classGroup1 = NULL,
                     classGroup2 = NULL, analysisName, groupName1,
                     groupName2, covariates = NULL, onlyPos = FALSE,
                     log2fcThreshold = 0.25, fdrThreshold = 0.05,
                     overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, index1, index2, class,
                                 classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)

    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    subsetIdx <- (ix1 | ix2)
    conditions <- rep(NA, ncol(inSCE))
    conditions[ix1] <- 'cond1'
    conditions[ix2] <- 'cond2'
    conditions <- conditions[!is.na(conditions)]
    annotData <- data.frame(condition = conditions,
                            row.names = colnames(inSCE)[subsetIdx])
    cov <- SummarizedExperiment::colData(inSCE)[subsetIdx, covariates,
                                                drop = FALSE]
    annotData <- cbind(annotData, cov)
    if (is.null(covariates)) {
      mod <- stats::model.matrix(~condition, annotData)
      mod0 <- stats::model.matrix(~1, annotData)
    } else {
      mod <- stats::model.matrix(
        stats::as.formula(paste0("~", paste0(c("condition", covariates),
                                             collapse = "+"))),
        data = annotData)
      mod0 <- stats::model.matrix(
        stats::as.formula(paste0("~", paste0(covariates, collapse = "+"))),
        data = annotData)
    }
    dat <- as.matrix(SummarizedExperiment::assay(inSCE[,subsetIdx], useAssay))
    if(!inherits(dat, 'matrix')){
        dat <- as.matrix(dat)
    }
    dat <- featureNameDedup(dat)
    n <- dim(dat)[2]
    m <- dim(dat)[1]
    df1 <- dim(mod)[2]
    df0 <- dim(mod0)[2]
    p <- rep(0, m)
    Id <- diag(n)
    resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
                        t(mod))
    resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*%
                         t(mod0))
    rss1 <- resid ^ 2 %*% rep(1, n)
    rss0 <- resid0 ^ 2 %*% rep(1, n)
    fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
    p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
    deg <- data.frame(row.names = rownames(dat), p.value = p,
                      padj = stats::p.adjust(p, method = 'fdr'))
    deg <- cbind(data.frame(Gene = rownames(deg),
                            Log2_FC = NA, stringsAsFactors = FALSE),
                 deg)
    rownames(deg) <- NULL
    colnames(deg) <- c("Gene", "Log2_FC", "Pvalue", "FDR")
    # Result Filtration
    if(isTRUE(onlyPos)){
      warning("ANOVA method does not generate log2FC value. `onlyPos` ",
              "argument ignored")
      #deg <- deg[deg$Log2_FC > 0,]
    }
    if(!is.null(fdrThreshold)){
      deg <- deg[deg$FDR < fdrThreshold,]
    }
    if(!is.null(log2fcThreshold)){
      warning("ANOVA method does not generate log2FC value. `log2fcThreshold` ",
              "argument ignored")
      #deg <- deg[abs(deg$Log2_FC) > log2fcThreshold,]
    }
    # Format output
    deg <- deg[order(deg$FDR, na.last = TRUE),]
    resultList$result <- deg
    resultList$method <- 'ANOVA'

    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
      S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    } else {
      S4Vectors::metadata(inSCE)$diffExp <- list()
      S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    }
    return(inSCE)
}

#' Perform differential expression analysis on SCE with MAST
#'
#' Condition specification allows two methods:
#' 1. Index level selection. Arguments \code{index1} and \code{index2} will be
#' used.
#' 2. Annotation level selection. Arguments \code{class}, \code{classGroup1} and
#' \code{classGroup2} will be used.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string specifying which assay to use for MAST
#' Default \code{"logcounts"}.
#' @param index1 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. Specifies
#' which cells are of interests. Default \code{NULL}.
#' @param index2 Any type of indices that can subset a
#' \linkS4class{SingleCellExperiment} inherited object by cells. specifies
#' the control group against those specified by \code{index1}. If
#' \code{NULL} when using index specification, \code{index1} cells will be
#' compared with all other cells. Default \code{NULL}.
#' @param class A vector/factor with \code{ncol(inSCE)} elements, or a character
#' scalar that specifies a column name of \code{colData(inSCE)}. Default
#' \code{NULL}.
#' @param classGroup1 a vector specifying which "levels" given in \code{class}
#' are of interests. Default \code{NULL}.
#' @param classGroup2 a vector specifying which "levels" given in \code{class}
#' is the control group against those specified by \code{classGroup1}. If
#' \code{NULL} when using annotation specification, \code{classGroup1} cells
#' will be compared with all other cells.
#' @param analysisName A character scalar naming the DEG analysis. Required
#' @param groupName1 A character scalar naming the group of interests. Required.
#' @param groupName2 A character scalar naming the control group. Required.
#' @param covariates A character vector of additional covariates to use when
#' building the model. All covariates must exist in
#' \code{names(colData(inSCE))}. Default \code{NULL}.
#' @param onlyPos Whether to only output DEG with positive log2_FC value.
#' Default \code{FALSE}.
#' @param log2fcThreshold Only out put DEGs with the absolute values of log2FC
#' greater than this value. Default \code{0.25}
#' @param fdrThreshold Only out put DEGs with FDR value less than this
#' value. Default \code{0.05}
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @param check_sanity Logical scalar. Whether to perform MAST's sanity check
#' to see if the counts are logged. Default \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' logcounts(sce) <- log(counts(sce) + 1)
#' sce <- runMAST(inSCE = sce, groupName1 = "Sample1",
#'  groupName2 = "Sample2", index1 = 1:20, index2 = 21:40,
#'  analysisName = "MAST")
#'
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$diffExp} updated with the results: a list named by
#' \code{analysisName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} storing the assay name that was used for
#' calculation, \code{$select} storing the cell selection indices (logical) for
#' each condition, \code{$result} storing a \code{\link{data.frame}} of
#' the DEGs summary, and \code{$method} storing \code{"MAST"}.
#' @export
runMAST <- function(inSCE, useAssay = 'logcounts', index1 = NULL,
                    index2 = NULL, class = NULL, classGroup1 = NULL,
                    classGroup2 = NULL, analysisName, groupName1,
                    groupName2, covariates = NULL, onlyPos = FALSE,
                    log2fcThreshold = NULL, fdrThreshold = 0.05,
                    overwrite = FALSE, check_sanity = TRUE){
    resultList <- .formatDEAList(inSCE, useAssay, index1, index2, class,
                                 classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)

    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    cells1 <- which(ix1)
    cells2 <- which(ix2)
    subsetIdx <- (ix1 | ix2)
    ## Extract
    mat <-
      SummarizedExperiment::assay(inSCE[,subsetIdx], useAssay)
    if(!inherits(mat, 'matrix')){
      mat <- as.matrix(mat)
    }
    mat <- featureNameDedup(mat)
    cond <- rep(NA, ncol(inSCE))
    cond[ix1] <- 'c1'
    cond[ix2] <- 'c2'
    cond <- cond[subsetIdx]
    cdat <- data.frame(wellKey = colnames(mat),
                       condition = as.factor(cond),
                       ngeneson = rep("", (length(cells1) + length(cells2))),
                       stringsAsFactors = FALSE)
    covariateDat <-
      data.frame(SummarizedExperiment::colData(inSCE[,subsetIdx])[,
                                                                  covariates,
                                                                  drop = FALSE])
    cdat <- cbind(cdat, covariateDat)
    sca <- MAST::FromMatrix(mat, cdat, check_sanity = check_sanity)
    cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
    SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
    cond <- factor(SummarizedExperiment::colData(sca)$condition)
    cond <- stats::relevel(cond, "c2")
    SummarizedExperiment::colData(sca)$condition <- cond
    # Calculation
    sca <- sca[which(MAST::freq(sca) > 0),]
    invisible(utils::capture.output(thresh <-
        MAST::thresholdSCRNACountMatrix(SummarizedExperiment::assay(sca),
                                        nbins = 20, min_per_bin = 30)))
    SummarizedExperiment::assays(sca) <-
        list(thresh = thresh$counts_threshold,
             tpm = SummarizedExperiment::assay(sca))
    SummarizedExperiment::colData(sca)$cngeneson <-
        scale(colSums(SummarizedExperiment::assay(sca) > 0))
    if(all(is.na(SummarizedExperiment::colData(sca)$cngeneson))){
        SummarizedExperiment::colData(sca)$cngeneson <- 0
    }
    zlmCond <- MAST::zlm(
        stats::as.formula(
            paste0("~", paste(c("condition", "cngeneson", covariates),
                              collapse = '+'))),
        sca)
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
    fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$FDR, na.last = TRUE),
    ]
    # Format output
    resultList$result <- as.data.frame(fcHurdleSig)
    resultList$method <- 'MAST'

    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    } else {
        S4Vectors::metadata(inSCE)$diffExp <- list()
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    }
    return(inSCE)
}
