#' Helper function for differential expression analysis methods that accepts
#' multiple ways of conditional subsetting and returns stable index format.
#' Meanwhile it does all the input checkings.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object. Required.
#' @param useAssay character. A string specifying which assay to use. Required.
#' @param useReducedDim character. A string specifying which reducedDim to use
#' for DE analysis. Usually a pathway analysis result matrix. Set 
#' \code{useAssay} to \code{NULL} when using. Required. 
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
.formatDEAList <- function(inSCE, useAssay, useReducedDim, index1 = NULL, index2 = NULL,
                           class = NULL, classGroup1 = NULL,
                           classGroup2 = NULL, groupName1, groupName2,
                           analysisName, covariates = NULL, overwrite = FALSE){
    # Input checks
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('"inSCE" should be a SingleCellExperiment inherited Object.')
    }
    if (is.null(useAssay) & is.null(useReducedDim)) {
      stop("`useAssay` or `useReducedDim` must be specified.")
    } else if (!is.null(useAssay) & !is.null(useReducedDim)) {
      stop("Only one of `useAssay` or `useReducedDim` can be specified.")
    } else {
      if (!is.null(useAssay)) {
        if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
          stop(paste('"useAssay" name: ', useAssay, ' not found.'))
        }
      } else {
        if(!useReducedDim %in% SingleCellExperiment::reducedDimNames(inSCE)){
          stop(paste('"useReducedDim" name: ', useReducedDim, ' not found.'))
        }
      }
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
            if(!isTRUE(overwrite)){
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
    isec <- intersect(which(ix1), which(ix2))
    if (length(isec) > 0) {
      stop("Cell(s) selected for both conditions: ",
           paste(colnames(inSCE)[isec], collapse = ", "))
    }
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
             useReducedDim = useReducedDim,
             groupNames = groupNames,
             select = select,
             annotation = annotation)
    )
}

# Filter the formated DEG table
# Could be used in runDE, getTable and plotDE functions
.filterDETable <- function(deg, onlyPos = FALSE, log2fcThreshold = NULL, 
                           fdrThreshold = NULL, minGroup1MeanExp = NULL,
                           maxGroup2MeanExp = NULL, minGroup1ExprPerc = NULL,
                           maxGroup2ExprPerc = NULL) {
  if (isTRUE(onlyPos)) {
    deg <- deg[deg$Log2_FC > 0,]
  }
  if (!is.null(fdrThreshold)) {
    deg <- deg[deg$FDR < fdrThreshold,]
  }
  if (!is.null(log2fcThreshold)) {
    deg <- deg[abs(deg$Log2_FC) >= log2fcThreshold,]
  }
  if (!is.null(minGroup1MeanExp)) {
    deg <- deg[deg$group1MeanExp >= minGroup1MeanExp,]
  }
  if (!is.null(maxGroup2MeanExp)) {
    deg <- deg[deg$group2MeanExp <= maxGroup2MeanExp,]
  }
  if (!is.null(minGroup1ExprPerc)) {
    deg <- deg[deg$group1ExprPerc >= minGroup1ExprPerc,]
  }
  if (!is.null(maxGroup2ExprPerc)) {
    deg <- deg[deg$group2ExprPerc <= maxGroup2ExprPerc,]
  }
  # Format output
  deg <- deg[order(deg$FDR, na.last = TRUE),]
  if (length(which(rowSums(is.na(deg)) > 2)) > 0) {
    deg <- deg[-which(rowSums(is.na(deg)) > 2),]
  }
  return(deg)
}

# Calculate extra metrics for each deg
# Including the mean expression value in group1,
# Percentage of cells in group1 that express each gene
# Percentage of cells in group2 that express each gene
.calculateDEMetrics <- function(deg, mat, ix1, ix2) {
  geneIdx <- rownames(mat) %in% deg$Gene
  # Mean expression in group1
  meanExp1 <- rowMeans(mat[geneIdx, ix1])
  meanExp1 <- meanExp1[deg$Gene]
  deg$group1MeanExp <- meanExp1
  # Mean expression in group2
  meanExp2 <- rowMeans(mat[geneIdx, ix2])
  meanExp2 <- meanExp2[deg$Gene]
  deg$group2MeanExp <- meanExp2
  # Expressed percentage in group1
  group1ExprPerc <- rowMeans(mat[geneIdx, ix1] > 0)
  group1ExprPerc <- group1ExprPerc[deg$Gene]
  deg$group1ExprPerc <-group1ExprPerc
  # Expressed percentage in group2
  group2ExprPerc <- rowMeans(mat[geneIdx, ix2] > 0)
  group2ExprPerc <- group2ExprPerc[deg$Gene]
  deg$group2ExprPerc <-group2ExprPerc
  return(deg)
}

#' Perform differential expression analysis on SCE object
#' @rdname runDEAnalysis
#' @details 
#' SCTK provides Limma, MAST, DESeq2, ANOVA and Wilcoxon test for differential
#' expression analysis, where DESeq2 expects non-negtive integer assay input 
#' while others expect logcounts. 
#' 
#' Condition specification allows two methods:
#' 1. Index level selection. Arguments \code{index1} and \code{index2} will be
#' used.
#' 2. Annotation level selection. Arguments \code{class}, \code{classGroup1} and
#' \code{classGroup2} will be used.
#' @seealso See \code{\link{plotDEGHeatmap}}, \code{\link{plotDEGRegression}}, 
#' \code{\link{plotDEGViolin}} and \code{\link{plotDEGVolcano}} for 
#' visualization method after running DE analysis. 
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param method Character. Specify which method to use when using 
#' \code{runDEAnalysis()}. Choose from \code{"wilcox"}, \code{"MAST"}, 
#' \code{"DESeq2"}, \code{"Limma"}, \code{"ANOVA"}. Default \code{"wilcox"}.
#' @param useAssay character. A string specifying which assay to use for the
#' DE regression. Default \code{"counts"} for DESeq2, \code{"logcounts"} for 
#' other methods. 
#' @param useReducedDim character. A string specifying which reducedDim to use
#' for DE analysis. Usually a pathway analysis result matrix. Set 
#' \code{useAssay} to \code{NULL} when using. Default \code{NULL}. 
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
#' greater than this value. Default \code{NULL}.
#' @param fdrThreshold Only out put DEGs with FDR value less than this
#' value. Default \code{NULL}.
#' @param minGroup1MeanExp Only out put DEGs with mean expression in group1 
#' greater then this value. Default \code{NULL}.
#' @param maxGroup2MeanExp Only out put DEGs with mean expression in group2 
#' less then this value. Default \code{NULL}.
#' @param minGroup1ExprPerc Only out put DEGs expressed in greater then this 
#' fraction of cells in group1. Default \code{NULL}.
#' @param maxGroup2ExprPerc Only out put DEGs expressed in less then this 
#' fraction of cells in group2. Default \code{NULL}.
#' @param overwrite A logical scalar. Whether to overwrite result if exists.
#' Default \code{FALSE}.
#' @param fullReduced Logical, DESeq2 only argument. Whether to apply LRT 
#' (Likelihood ratio test) with a 'full' model. Default \code{TRUE}.
#' @param check_sanity Logical, MAST only argument. Whether to perform MAST's 
#' sanity check to see if the counts are logged. Default \code{TRUE}.
#' @param ... Arguments to pass to specific methods when using the generic 
#' \code{runDEAnalysis()}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' sce <- runDEAnalysis(method = "Limma", inSCE = sce, groupName1 = "group1",
#'  groupName2 = "group2", index1 = seq(20), index2 = seq(21,40),
#'  analysisName = "Limma")
#'
#' @return The input \linkS4class{SingleCellExperiment} object with
#' \code{metadata(inSCE)$diffExp} updated with the results: a list named by
#' \code{analysisName}, with \code{$groupNames} containing the naming of the
#' two conditions, \code{$useAssay} and \code{$useReducedDim} storing the matrix
#' name that was used for calculation, \code{$select} storing the cell selection
#' indices (logical) for each condition, \code{$result} storing a 
#' \code{\link{data.frame}} of the DEGs summary, and \code{$method} storing the 
#' character method name used.
#' @export
runDEAnalysis <- function(method = c('wilcox', 'MAST', 'DESeq2', 'Limma', 
                                     'ANOVA'), ...){
    method <- match.arg(method)
    funcList <- list(MAST = runMAST,
                     DESeq2 = runDESeq2,
                     Limma = runLimmaDE,
                     ANOVA = runANOVA,
                     wilcox = runWilcox)
    funcList[[method]](...)
}

#' @rdname runDEAnalysis
#' @export
runDESeq2 <- function(inSCE, useAssay = 'counts', useReducedDim = NULL, 
                      index1 = NULL, index2 = NULL, class = NULL, 
                      classGroup1 = NULL, classGroup2 = NULL, analysisName, 
                      groupName1, groupName2, covariates = NULL, 
                      fullReduced = TRUE, onlyPos = FALSE, 
                      log2fcThreshold = NULL, fdrThreshold = NULL, 
                      minGroup1MeanExp = NULL, maxGroup2MeanExp = NULL, 
                      minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL,
                      overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, useReducedDim, index1, index2, 
                                 class, classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)
    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    subsetIdx <- (ix1 | ix2)
    conditions <- rep(NA, ncol(inSCE))
    conditions[ix1] <- 'cond1'
    conditions[ix2] <- 'cond2'
    conditions <- conditions[!is.na(conditions)]
    annotData <- data.frame(condition = factor(conditions),
                            row.names = colnames(inSCE)[subsetIdx])
    cov <- SummarizedExperiment::colData(inSCE)[subsetIdx, covariates,
                                                drop = FALSE]
    annotData <- cbind(annotData, cov)
    if (!is.null(useAssay)) {
      mat <- expData(inSCE[,subsetIdx], useAssay)
    } else {
      mat <- t(expData(inSCE[,subsetIdx], useReducedDim))
    }
    
    if(!inherits(mat, 'matrix')){
        mat <- as.matrix(mat)
    }
    if (any(duplicated(rownames(mat)))) {
      warning("Duplicated feature names found in given dataset. Making them ",
              "unique in the result. They will not show in plots.")
      mat <- dedupRowNames(mat)
    }
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
    deg <- cbind(data.frame(Gene = as.character(rownames(deg)),
                            stringsAsFactors = FALSE),
                 deg)
    rownames(deg) <- NULL
    colnames(deg) <- c('Gene', 'Log2_FC', 'Pvalue', 'FDR')
    deg$Log2_FC <- - deg$Log2_FC # JNMLP
    deg <- .calculateDEMetrics(deg, expData(inSCE, useAssay), ix1, ix2)
    deg <- .filterDETable(deg, onlyPos, log2fcThreshold, fdrThreshold, 
                          minGroup1MeanExp, maxGroup2MeanExp, 
                          minGroup1ExprPerc, maxGroup2ExprPerc)
    
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

#' @rdname runDEAnalysis
#' @export
runLimmaDE <- function(inSCE, useAssay = 'logcounts', useReducedDim = NULL, 
                       index1 = NULL, index2 = NULL, class = NULL, 
                       classGroup1 = NULL, classGroup2 = NULL, analysisName, 
                       groupName1, groupName2, covariates = NULL, 
                       onlyPos = FALSE, log2fcThreshold = NULL, 
                       fdrThreshold = NULL, minGroup1MeanExp = NULL, 
                       maxGroup2MeanExp = NULL, minGroup1ExprPerc = NULL, 
                       maxGroup2ExprPerc = NULL, overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, useReducedDim, index1, index2,
                                 class, classGroup1, classGroup2, groupName1,
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
    if (!is.null(useAssay)) {
      mat <- expData(inSCE[,subsetIdx], useAssay)
    } else {
      mat <- t(expData(inSCE[,subsetIdx], useReducedDim))
    }
    
    if(!inherits(mat, 'matrix')){
        mat <- as.matrix(mat)
    }
    if (any(duplicated(rownames(mat)))) {
      warning("Duplicated feature names found in given dataset. Making them ",
              "unique in the result. They will not show in plots.")
      mat <- dedupRowNames(mat)
    }
    design <- stats::model.matrix(
        stats::as.formula(paste0("~", paste0(c('condition', covariates),
                                             collapse = "+"))),
        data = annotData)
    fit <- limma::lmFit(mat, design)
    ebayes <- limma::eBayes(fit)
    coef <- seq(ncol(ebayes))[-1]
    deg <- limma::topTable(ebayes, coef = coef, adjust = 'fdr',
                           number = nrow(inSCE))
    deg <- deg[,c(-2, -3, -6)]
    deg <- cbind(data.frame(Gene = as.character(rownames(deg)),
                            stringsAsFactors = FALSE),
                 deg)
    rownames(deg) <- NULL
    colnames(deg) <- c("Gene", "Log2_FC", "Pvalue", "FDR")
    deg$Log2_FC <- - deg$Log2_FC # JNMLP
    deg <- .calculateDEMetrics(deg, expData(inSCE, useAssay), ix1, ix2)
    deg <- .filterDETable(deg, onlyPos, log2fcThreshold, fdrThreshold, 
                          minGroup1MeanExp, maxGroup2MeanExp, 
                          minGroup1ExprPerc, maxGroup2ExprPerc)

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

#' @rdname runDEAnalysis
#' @export
runANOVA <- function(inSCE, useAssay = 'logcounts', useReducedDim = NULL, 
                     index1 = NULL, index2 = NULL, class = NULL, 
                     classGroup1 = NULL, classGroup2 = NULL, analysisName, 
                     groupName1, groupName2, covariates = NULL, onlyPos = FALSE,
                     log2fcThreshold = NULL, fdrThreshold = NULL, 
                     minGroup1MeanExp = NULL, maxGroup2MeanExp = NULL, 
                     minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL,
                     overwrite = FALSE){
    resultList <- .formatDEAList(inSCE, useAssay, useReducedDim, index1, index2,
                                 class, classGroup1, classGroup2, groupName1,
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
    if (!is.null(useAssay)) {
      dat <- expData(inSCE[,subsetIdx], useAssay)
    } else {
      dat <- t(expData(inSCE[,subsetIdx], useReducedDim))
    }
    if(!inherits(dat, 'matrix')){
        dat <- as.matrix(dat)
    }
    if (any(duplicated(rownames(dat)))) {
      warning("Duplicated feature names found in given dataset. Making them ",
              "unique in the result. They will not show in plots.")
      dat <- dedupRowNames(dat)
    }

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
    deg <- data.frame(Gene = as.character(rownames(dat)),
                      Log2_FC =NA, Pvalue = p,
                      FDR = stats::p.adjust(p, method = 'fdr'),
                      stringsAsFactors = FALSE)
    rownames(deg) <- NULL
    
    if (!is.null(useAssay)) {
      cond1.assay <- expData(inSCE, useAssay)[, ix1]
      cond2.assay <- expData(inSCE, useAssay)[, ix2]
    } else {
      cond1.assay <- t(expData(inSCE, useReducedDim))[, ix1]
      cond2.assay <- t(expData(inSCE, useReducedDim))[, ix2]
    }
    # Assuming that useAssay is log-normalized counts
    deg$Log2_FC <- rowMeans(cond1.assay) - rowMeans(cond2.assay)
    deg <- .calculateDEMetrics(deg, expData(inSCE, useAssay), ix1, ix2)
    deg <- .filterDETable(deg, onlyPos, log2fcThreshold, fdrThreshold,
                          minGroup1MeanExp, maxGroup2MeanExp, 
                          minGroup1ExprPerc, maxGroup2ExprPerc)
    
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

#' @rdname runDEAnalysis
#' @export
runMAST <- function(inSCE, useAssay = 'logcounts', useReducedDim = NULL,
                    index1 = NULL, index2 = NULL, class = NULL, 
                    classGroup1 = NULL, classGroup2 = NULL, analysisName, 
                    groupName1, groupName2, covariates = NULL, onlyPos = FALSE,
                    log2fcThreshold = NULL, fdrThreshold = NULL,
                    minGroup1MeanExp = NULL, maxGroup2MeanExp = NULL, 
                    minGroup1ExprPerc = NULL, maxGroup2ExprPerc = NULL,
                    overwrite = FALSE, check_sanity = TRUE){
    resultList <- .formatDEAList(inSCE, useAssay, useReducedDim, index1, index2, 
                                 class, classGroup1, classGroup2, groupName1,
                                 groupName2, analysisName, covariates,
                                 overwrite)

    ix1 <- resultList$select$ix1
    ix2 <- resultList$select$ix2
    cells1 <- which(ix1)
    cells2 <- which(ix2)
    subsetIdx <- (ix1 | ix2)
    ## Extract
    if (!is.null(useAssay)) {
      mat <- expData(inSCE[,subsetIdx], useAssay)
    } else {
      mat <- t(expData(inSCE[,subsetIdx], useReducedDim))
    }
    
    if(!inherits(mat, 'matrix')){
      mat <- as.matrix(mat)
    }
    if (any(duplicated(rownames(mat)))) {
      warning("Duplicated feature names found in given dataset. Making them ",
              "unique in the result. They will not show in plots.")
      mat <- dedupRowNames(mat)
    }
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
    fcHurdleSig <- fcHurdleSig[, -c(4, 5)]
    names(fcHurdleSig)[c(1, 2, 3, 4)] <- c("Gene", "Pvalue",
                                           "Log2_FC", "FDR")
    fcHurdleSig$Gene <- as.character(fcHurdleSig$Gene)
    fcHurdleSig <- .calculateDEMetrics(fcHurdleSig, expData(inSCE, useAssay), 
                                       ix1, ix2)
    fcHurdleSig <- .filterDETable(fcHurdleSig, onlyPos, log2fcThreshold, 
                                  fdrThreshold, minGroup1MeanExp, 
                                  maxGroup2MeanExp, minGroup1ExprPerc, 
                                  maxGroup2ExprPerc)

    resultList$result <- fcHurdleSig
    resultList$method <- 'MAST'

    if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    } else {
        S4Vectors::metadata(inSCE)$diffExp <- list()
        S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
    }
    return(inSCE)
}

#' @rdname runDEAnalysis
#' @export
runWilcox <- function(inSCE, useAssay = 'logcounts', useReducedDim = NULL, 
                      index1 = NULL, index2 = NULL, class = NULL, 
                      classGroup1 = NULL, classGroup2 = NULL, analysisName, 
                      groupName1, groupName2, covariates = NULL, 
                      onlyPos = FALSE, log2fcThreshold = NULL, 
                      fdrThreshold = NULL, minGroup1MeanExp = NULL, 
                      maxGroup2MeanExp = NULL, minGroup1ExprPerc = NULL, 
                      maxGroup2ExprPerc = NULL,overwrite = FALSE){
  resultList <- .formatDEAList(inSCE, useAssay, useReducedDim, index1, index2, 
                               class, classGroup1, classGroup2, groupName1,
                               groupName2, analysisName, covariates,
                               overwrite)
  if (!is.null(covariates)) {
    warning("Wilcoxon test from Scran does not support covariate modeling.
            Ignoring this argument in calculation, but will be included in
            plotting.")
  }
  ix1 <- resultList$select$ix1
  ix2 <- resultList$select$ix2
  conditions <- rep(NA, ncol(inSCE))
  conditions[ix1] <- 'cond1'
  conditions[ix2] <- 'cond2'
  if (!is.null(useAssay)) {
    mat <- expData(inSCE, useAssay)
  } else {
    mat <- t(expData(inSCE, useReducedDim))
  }
  result <- scran::pairwiseWilcox(mat, groups = conditions)
  # result <- scran::pairwiseWilcox(inSCE, groups = conditions,
  #                                 assay.type = useAssay)
  if (result$pairs$first[1] == "cond1") {
    table <- result$statistics[[1]]
  } else {
    table <- result$statistics[[2]]
  }
  # Generate LogFC value
  
  if (!is.null(useAssay)) {
    cond1.assay <- expData(inSCE, useAssay)[rownames(table), ix1]
    cond2.assay <- expData(inSCE, useAssay)[rownames(table), ix2]
  } else {
    cond1.assay <- t(expData(inSCE, useReducedDim))[rownames(table), ix1]
    cond2.assay <- t(expData(inSCE, useReducedDim))[rownames(table), ix2]
  }
  # Assuming that useAssay is log-normalized counts
  table$Log2_FC <- rowMeans(cond1.assay) - rowMeans(cond2.assay)
  
  deg <- data.frame(Gene = as.character(rownames(table)),
                    Log2_FC = table$Log2_FC,
                    Pvalue = table$p.value,
                    FDR = table$FDR)
  rownames(deg) <- NULL
  deg <- .calculateDEMetrics(deg, mat, ix1, ix2)
  deg <- .filterDETable(deg, onlyPos, log2fcThreshold, fdrThreshold,
                        minGroup1MeanExp, maxGroup2MeanExp, 
                        minGroup1ExprPerc, maxGroup2ExprPerc)
  
  resultList$result <- deg
  resultList$method <- 'wilcox'
  if ("diffExp" %in% names(S4Vectors::metadata(inSCE))){
    S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
  } else {
    S4Vectors::metadata(inSCE)$diffExp <- list()
    S4Vectors::metadata(inSCE)$diffExp[[analysisName]] <- resultList
  }
  return(inSCE)
}
