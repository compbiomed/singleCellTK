#' Estimate numbers of detected genes, significantly differentially expressed
#' genes, and median significant effect size
#'
#' @param originalData \linkS4class{SingleCellExperiment} object storing all
#' assay data from the shiny app.
#' @param useAssay Character. The name of the assay to be used for subsampling.
#' @param minCount Numeric. The minimum number of reads found for a gene to be
#' considered detected.
#' @param minCells Numeric. The minimum number of cells a gene must have at
#' least 1 read in for it to be considered detected.
#' @param maxDepth Numeric. The highest number of total reads to be simulated.
#' @param realLabels Character. The name of the condition of interest. Must match
#' a name from sample data.
#' @param depthResolution Numeric. How many different read depth should the
#' script simulate? Will simulate a number of experimental designs ranging from
#' 10 reads to maxReadDepth, with logarithmic spacing.
#' @param iterations Numeric. How many times should each experimental design be
#' simulated?
#'
#' @return A 3-dimensional array, with dimensions = c(iterations,
#' depthResolution, 3). [,,1] contains the number of detected genes in each
#' simulated dataset, [,,2] contains the number of significantly differentially
#' expressed genes in each simulation, and [,,3] contains the mediansignificant
#' effect size in each simulation. If no genes are significantly differentially
#' expressed, the median effect size defaults to infinity.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' subset <- mouseBrainSubsetSCE[1:1000,]
#' res <- downSampleDepth(subset,
#'                        realLabels = "level1class",
#'                        iterations=2)
#'
downSampleDepth <- function(originalData, useAssay = "counts", minCount = 10, minCells = 3,
                            maxDepth = 10000000, realLabels,
                            depthResolution = 10, iterations = 10){
  realLabels <- SingleCellExperiment::colData(originalData)[, realLabels]
  originalData <- SummarizedExperiment::assay(originalData, useAssay)
  foundGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  minEffectSizeMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  numSigGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  depths <- floor(10 ^ (seq(0, log10(maxDepth), length.out = depthResolution)))
  cells <- dim(originalData)[2]
  effectSizes <- calcEffectSizes(originalData, realLabels)
  for (i in seq_len(depthResolution)){
    for (j in seq_len(iterations)){
      tempData <- generateSimulatedData(totalReads = depths[i], cells,
                                        as.matrix(originalData),
                                        realLabels = as.factor(realLabels))
      tempSigDiff <- subDiffEx(tempData)
      foundGenesMatrix[j, i] <- sum(apply(tempData[-1, ], 1, function(x){
        sum(x > 0) >= minCells && sum(x) >= minCount
      }))
      numSigGenesMatrix[j, i] <- sum(tempSigDiff <= 0.05)
      minEffectSizeMatrix[j, i] <-
        abs(min(abs(effectSizes[which(tempSigDiff <= 0.05)])))
    }
  }
  outArray <- array(c(foundGenesMatrix, minEffectSizeMatrix, numSigGenesMatrix),
                    dim = c(iterations, depthResolution, 3))
  return(outArray)
}

#' Estimate numbers of detected genes, significantly differentially expressed
#' genes, and median significant effect size
#'
#' @param originalData The \linkS4class{SingleCellExperiment} object storing all
#' assay data from the shiny app.
#' @param useAssay Character. The name of the assay to be used for subsampling.
#' @param minCountDetec Numeric. The minimum number of reads found for a gene to
#' be considered detected.
#' @param minCellsDetec Numeric. The minimum number of cells a gene must have at
#' least 1 read in for it to be considered detected.
#' @param minCellnum Numeric. The minimum number of virtual cells to include in
#' the smallest simulated dataset.
#' @param maxCellnum Numeric. The maximum number of virtual cells to include in
#' the largest simulated dataset
#' @param realLabels Character. The name of the condition of interest. Must match
#' a name from sample data. If only two factors present in the corresponding
#' colData, will default to t-test. If multiple factors, will default to ANOVA.
#' @param depthResolution Numeric. How many different read depth should the
#' script simulate? Will simulate a number of experimental designs ranging from
#' 10 reads to maxReadDepth, with logarithmic spacing.
#' @param iterations Numeric. How many times should each experimental design be
#' simulated?
#' @param totalReads Numeric. How many aligned reads to put in each simulated
#' dataset.
#'
#' @return A 3-dimensional array, with dimensions = c(iterations,
#' depthResolution, 3). [,,1] contains the number of detected genes in each
#' simulated dataset, [,,2] contains the number of significantly differentially
#' expressed genes in each simulation, and [,,3] contains the mediansignificant
#' effect size in each simulation. If no genes are significantly differentially
#' expressed, the median effect size defaults to infinity.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' subset <- mouseBrainSubsetSCE[1:1000,]
#' res <- downSampleCells(subset,
#'                        realLabels = "level1class",
#'                        iterations=2)
downSampleCells <- function(originalData, useAssay = "counts", minCountDetec = 10, minCellsDetec = 3,
                            minCellnum = 10, maxCellnum = 1000, realLabels,
                            depthResolution = 10, iterations = 10,
                            totalReads = 1000000){
  realLabels <- SingleCellExperiment::colData(originalData)[, realLabels]
  originalData <- SummarizedExperiment::assay(originalData, useAssay)
  foundGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  minEffectSizeMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  numSigGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  cellNums <- seq(minCellnum, maxCellnum, length.out = depthResolution)
  effectSizes <- calcEffectSizes(originalData, realLabels)
  for (i in seq_len(depthResolution)){
    for (j in seq_len(iterations)){
      tempData <- generateSimulatedData(totalReads = totalReads,
                                        cells = cellNums[i],
                                        as.matrix(originalData),
                                        realLabels = as.factor(realLabels))
      tempSigDiff <- subDiffEx(tempData)
      foundGenesMatrix[j, i] <- sum(apply(tempData[-1, ], 1, function(x){
        sum(x > 0) >= minCellsDetec && sum(x) >= minCountDetec
      }))
      numSigGenesMatrix[j, i] <- sum(tempSigDiff <= 0.05)
      minEffectSizeMatrix[j, i] <-
        abs(min(abs(effectSizes[which(tempSigDiff <= 0.05)])))
    }
  }
  outArray <- array(c(foundGenesMatrix, minEffectSizeMatrix, numSigGenesMatrix),
                    dim = c(iterations, depthResolution, 3))
  return(outArray)
}

#' Generates a single simulated dataset, bootstrapping from the input counts
#' matrix.
#'
#' @param totalReads Numeric. The total number of reads in the simulated
#' dataset, to be split between all simulated cells.
#' @param originalData Matrix. The original raw read count matrix. When used
#' within the Shiny app, this will be assay(SCEsetObject, "counts").
#' @param cells Numeric. The number of virtual cells to simulate.
#' @param realLabels Factor. The condition labels for differential expression.
#' If only two factors present, will default to t-test. If multiple factors,
#' will default to ANOVA.
#'
#' @return A simulated counts matrix, the first row of which contains the 'true'
#' labels for each virtual cell.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- generateSimulatedData(
#'          totalReads = 1000, cells=10,
#'          originalData = assay(mouseBrainSubsetSCE, "counts"),
#'          realLabels = colData(mouseBrainSubsetSCE)[, "level1class"])
#'
generateSimulatedData <- function(totalReads, cells, originalData, realLabels){
  cells <- sample(dim(originalData)[2], size = cells, replace = TRUE)
  totalReads <- floor(totalReads / length(cells))
  originalData <- t(t(originalData) / apply(originalData, 2, sum))
  output <- matrix(nrow = dim(originalData)[1], ncol = length(cells))
  for (i in seq_along(cells)){
    output[, i] <- stats::rmultinom(1, totalReads, originalData[, cells[i]])
  }
  return(rbind(realLabels[cells], output))
}

#' Returns significance data from a snapshot.
#'
#' @param originalData The \linkS4class{SingleCellExperiment} object storing all
#' assay data from the shiny app.
#' @param useAssay Character. The name of the assay to be used for subsampling.
#' @param realLabels Character. The name of the condition of interest. Must match
#' a name from sample data.
#' @param totalReads Numeric. The total number of reads in the simulated
#' dataset, to be split between all simulated cells.
#' @param cells Numeric. The number of virtual cells to simulate.
#' @param iterations Numeric. How many times should each experimental design be
#' simulated.
#'
#' @return A matrix of significance information from a snapshot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- iterateSimulations(mouseBrainSubsetSCE, realLabels = "level1class",
#'                           totalReads = 1000, cells = 10, iterations = 2)
#'
iterateSimulations <- function(originalData, useAssay = "counts", realLabels, totalReads, cells,
                               iterations){
  realLabels <- SingleCellExperiment::colData(originalData)[, realLabels]
  originalData <- SummarizedExperiment::assay(originalData, useAssay)
  sigMatrix <- matrix(nrow = dim(originalData)[1])
  for (i in seq_len(iterations)){
    tempData <- generateSimulatedData(totalReads, cells, originalData,
                                      realLabels = as.factor(realLabels))
    sigMatrix <- cbind(sigMatrix, subDiffEx(tempData))
  }
  return(sigMatrix[, -1])
}

#' Passes the output of generateSimulatedData() to differential expression
#' tests, picking either t-tests or ANOVA for data with only two conditions or
#' multiple conditions, respectively.
#'
#' @param tempData Matrix. The output of generateSimulatedData(), where the
#' first row contains condition labels.
#'
#' @return subDiffEx(): A vector of fdr-adjusted p-values for all genes.
#' Nonviable results (such as for genes with 0 counts in a simulated dataset)
#' are coerced to 1.
#'
#' @describeIn subDiffEx 
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- generateSimulatedData(
#'          totalReads = 1000, cells=10,
#'          originalData = assay(mouseBrainSubsetSCE, "counts"),
#'          realLabels = colData(mouseBrainSubsetSCE)[, "level1class"])
#' tempSigDiff <- subDiffEx(res)
#'
subDiffEx <- function(tempData){
  realLabels <- tempData[1, ]
  output <- tempData[-1, ]
  if (length(levels(as.factor(realLabels))) > 2){
    fdr <- subDiffExANOVA(output, realLabels)
  } else if (length(levels(as.factor(realLabels))) == 2){
    fdr <- subDiffExttest(output, realLabels)
  }
  else{
    stop("Only 1 (or 0?) factor in ", levels(as.factor(realLabels)))
  }
  fdr[which(is.na(fdr))] <- 1
  return(fdr)
}

#' @describeIn subDiffEx Runs t-tests on all genes in a simulated dataset with 2
#' conditions, and adjusts for FDR.
#'
#' @param countMatrix Matrix. A simulated counts matrix, sans labels.
#' @param class.labels Factor. The condition labels for the simulated cells.
#' Will be coerced into 1's and 0's.
#' @param test.type Type of test to perform. The default is t.equalvar.
#'
#' @return subDiffExttest(): A vector of fdr-adjusted p-values for all genes.
#' Nonviable results (such as for genes with 0 counts in a simulated dataset)
#' are coerced to 1.
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #sort first 100 expressed genes
#' ord <- rownames(mouseBrainSubsetSCE)[
#'   order(rowSums(assay(mouseBrainSubsetSCE, "counts")),
#'         decreasing = TRUE)][1:100]
#' #subset to those first 100 genes
#' subset <- mouseBrainSubsetSCE[ord, ]
#' res <- generateSimulatedData(totalReads = 1000, cells=10,
#'                              originalData = assay(subset, "counts"),
#'                              realLabels = colData(subset)[, "level1class"])
#' realLabels <- res[1, ]
#' output <- res[-1, ]
#' fdr <- subDiffExttest(output, realLabels)
#'
subDiffExttest <- function(countMatrix, class.labels, test.type = "t.equalvar") {
  class.labels <- as.numeric(as.factor(class.labels))
  class.labels <- class.labels - 1
  class.labels[class.labels > 0] <- 1
  tval <- multtest::mt.teststat(countMatrix, classlabel = class.labels,
                                test = test.type, nonpara = "n")
  df <- (ncol(countMatrix) - 2)
  pval <- 2 * (1 - stats::pt(abs(tval), df))
  fdr <- stats::p.adjust(pval, method = "fdr")
  return(fdr)
}

#' @describeIn subDiffEx Runs ANOVA on all genes in a simulated dataset with
#' more than 2 conditions, and adjusts for FDR.
#'
#' @param condition Factor. The condition labels for the simulated cells.
#'
#' @return subDiffExANOVA(): A vector of fdr-adjusted p-values for all genes.
#' Nonviable results (such as for genes with 0 counts in a simulated dataset)
#' are coerced to 1.
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' #sort first 100 expressed genes
#' ord <- rownames(mouseBrainSubsetSCE)[
#'   order(rowSums(assay(mouseBrainSubsetSCE, "counts")),
#'         decreasing = TRUE)][1:100]
#' # subset to those first 100 genes
#' subset <- mouseBrainSubsetSCE[ord, ]
#' res <- generateSimulatedData(totalReads = 1000, cells=10,
#'                              originalData = assay(subset, "counts"),
#'                              realLabels = colData(subset)[, "level2class"])
#' realLabels <- res[1, ]
#' output <- res[-1, ]
#' fdr <- subDiffExANOVA(output, realLabels)
#'
subDiffExANOVA <- function(countMatrix, condition){
  mod <- stats::model.matrix(~as.factor(condition))
  mod0 <- stats::model.matrix(~1, data = as.factor(condition))
  dat <- log2(countMatrix + 1)
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
  fstats <- ((rss0 - rss1) / (df1 - df0)) / ((rss1 + 1e-50) / (n - df1))
  p <- 1 - stats::pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  return(stats::p.adjust(p, method = "fdr"))
}

#' Finds the effect sizes for all genes in the original dataset, regardless of
#' significance.
#'
#' @param countMatrix Matrix. A simulated counts matrix, sans labels.
#' @param condition Factor. The condition labels for the simulated cells. If
#' more than 2 conditions are given, the first will be compared to all others by
#' default.
#'
#' @return A vector of cohen's d effect sizes for each gene.
#'
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' res <- calcEffectSizes(assay(mouseBrainSubsetSCE, "counts"),
#'                        condition = colData(mouseBrainSubsetSCE)[, "level1class"])
#'
calcEffectSizes <- function(countMatrix, condition){
  groups <- levels(as.factor(unlist(condition)))
  return((apply(countMatrix[, condition == unlist(groups)[1]], 1, mean) -
            apply(countMatrix[, condition != unlist(groups)[1]], 1, mean)) /
           apply(countMatrix, 1, S4Vectors::sd))
}

powerCalc <- function(datamatrix, sampleSizeRange=c(1000, 10000000),
                      byBatch=FALSE, batch=NULL, numSize=25) {
  if (byBatch == FALSE){
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                numSize))
    for (j in seq_len(dim(datamatrix)[2])) {
      probs <- as.numeric(datamatrix[, j] / sum(datamatrix[, j]))
      for (i in seq_along(probs)){
        discoveryPower <- 1 - DelayedArray::dbinom(0, size = floor(seq.int(
          from = sampleSizeRange[1], to = sampleSizeRange[2],
          length.out = numSize)), prob = probs[i])
        outmat[i, j, ] <- discoveryPower
      }
    }
  }
  else {
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                numSize))
    for (j in seq_len(nlevels(batch))) {
      probs <- datamatrix[, which(batch == levels(batch)[j])] /
        sum(datamatrix[, which(batch == levels(batch)[j])])
      for (i in seq_along(as.vector(probs))) {
        discoveryPower <- 1 - DelayedArray::dbinom(0, size = floor(seq.int(
          from = sampleSizeRange[1], to = sampleSizeRange[2],
          length.out = numSize)), prob = as.vector(probs)[i])
        outmat[((i - 1) %% dim(datamatrix)[1]) + 1,
               which(batch == levels(batch)[j])[
                 ceiling(i / dim(datamatrix)[1])], ] <- discoveryPower
      }
    }
  }
  return(outmat)
}
