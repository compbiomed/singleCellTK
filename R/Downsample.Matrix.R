#' Estimate numbers of detected genes, significantly differentially expressed
#' genes, and median significant effect size
#'
#' @param originalData Matrix. The original raw readcount matrix. When used
#' within the Shiny app, this will be assay(SCEsetObject, "counts").
#' @param minCount Numeric. The minimum number of reads found for a gene to be
#' considered detected.
#' @param minCells Numeric. The minimum number of cells a gene must have at
#' least 1 read in for it to be considered detected.
#' @param maxDepth Numeric. The highest number of total reads to be simulated.
#' @param realLabels Factor. The condition labels for differential expression.
#' If only two factors present, will default to t-test. If multiple factors,
#' will default to ANOVA.
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
#'
DownsampleDepth <- function(originalData, minCount = 10, minCells = 3,
                            maxDepth = 10000000, realLabels,
                            depthResolution = 10, iterations = 10){
  realLabels <- colData(originalData)[, realLabels]
  originalData <- counts(originalData)
  foundGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  minEffectSizeMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  numSigGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  depths <- floor(10 ^ (seq(0, log10(maxDepth), length.out = depthResolution)))
  cells <- dim(originalData)[2]
  effectSizes <- calcEffectSizes(originalData, realLabels)
  for(i in 1:depthResolution){
    for(j in 1:iterations){
      tempData <- generateSimulatedData(totalReads = depths[i], cells,
                                        as.matrix(originalData),
                                        realLabels = as.factor(realLabels))
      tempSigDiff <- subDiffEx(tempData)
      foundGenesMatrix[j, i] <- sum(apply(tempData[-1, ], 1, function(x){
        sum(x > 0) >= minCells && sum(x) >= minCount
      }))
      numSigGenesMatrix[j, i] <- sum(tempSigDiff <= 0.05)
      minEffectSizeMatrix[j, i] <- abs(min(abs(effectSizes[which(tempSigDiff <= 0.05)])))
    }
  }
  outArray <- array(c(foundGenesMatrix, minEffectSizeMatrix, numSigGenesMatrix),
                    dim = c(iterations, depthResolution, 3))
  return(outArray)
}

#' Estimate numbers of detected genes, significantly differentially expressed
#' genes, and median significant effect size
#'
#' @param originalData Matrix. The original raw readcount matrix. When used
#' within the Shiny app, this will be assay(SCEsetObject, "counts").
#' @param minCountDetec Numeric. The minimum number of reads found for a gene to
#' be considered detected.
#' @param minCellsDetec Numeric. The minimum number of cells a gene must have at
#' least 1 read in for it to be considered detected.
#' @param minCellnum Numeric. The minimum number of virtual cells to include in
#' the smallest simulated dataset.
#' @param maxCellnum Numeric. The maximum number of virtual cells to include in
#' the largest simulated dataset
#' @param realLabels Factor. The condition labels for differential expression.
#' If only two factors present, will default to t-test. If multiple factors,
#' will default to ANOVA.
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
#'
DownsampleCells <- function(originalData, minCountDetec = 10, minCellsDetec = 3,
                            minCellnum = 10, maxCellnum = 1000, realLabels,
                            depthResolution = 10, iterations = 10,
                            totalReads = 1000000){
  realLabels <- colData(originalData)[, realLabels]
  originalData <- counts(originalData)
  foundGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  minEffectSizeMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  numSigGenesMatrix <- matrix(nrow = iterations, ncol = depthResolution)
  cellNums <- seq(minCellnum, maxCellnum, length.out = depthResolution)
  effectSizes <- calcEffectSizes(originalData, realLabels)
  for(i in 1:depthResolution){
    for(j in 1:iterations){
      tempData <- generateSimulatedData(totalReads = totalReads,
                                        cells = cellNums[i],
                                        as.matrix(originalData),
                                        realLabels = as.factor(realLabels))
      tempSigDiff <- subDiffEx(tempData)
      foundGenesMatrix[j, i] <- sum(apply(tempData[-1, ], 1, function(x){
        sum(x > 0) >= minCellsDetec && sum(x) >= minCountDetec
      }))
      numSigGenesMatrix[j, i] <- sum(tempSigDiff <= 0.05)
      minEffectSizeMatrix[j, i] <- abs(min(abs(effectSizes[which(tempSigDiff <= 0.05)])))
    }
  }
  outArray <- array(c(foundGenesMatrix, minEffectSizeMatrix, numSigGenesMatrix),
                    dim = c(iterations, depthResolution, 3))
  return(outArray)
}

#' Downsample Data
#'
#' @param datamatrix TODO:document
#' @param newcounts TODO:document
#' @param byBatch TODO:document
#' @param batch TODO:document
#' @param iterations TODO:document
#'
#' @return Downsampled matrix
#' @export
#'
Downsample <- function(datamatrix, newcounts = c(4, 16, 64, 256, 1024, 4096,
                                                 16384, 65536, 262144),
                       byBatch = FALSE, batch = NULL, iterations = 10) {
  if (byBatch == FALSE) {
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                length(newcounts), iterations))
    for (j in 1:dim(datamatrix)[2]) {
      probs <- datamatrix[, j] / sum(datamatrix[, j])
      for (k in 1:length(newcounts)) {
        samps <- stats::rmultinom(iterations, newcounts[k], probs)
        for (l in 1:iterations) {
          outmat[, j, k, l] <- samps[, l]
        }
      }
    }
  }
  else {
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                length(newcounts), iterations))
    for (j in 1:nlevels(batch)) {
      probs <- datamatrix[, which(batch == levels(batch)[j])] / sum(datamatrix[, which(batch == levels(batch)[j])])
      for (k in 1:length(newcounts)) {
        samps <- stats::rmultinom(iterations, newcounts[k], as.vector(probs))
        for (l in 1:iterations) {
          outmat[, which(batch == levels(batch)[j]), k, l] <- as.matrix(samps[, l], nrow = dim(datamatrix)[1])
        }
      }
    }
  }
  return(outmat)
}

#' Generates a single simulated dataset, bootstrapping from the input counts
#' matrix.
#'
#' @param originalData Matrix. The original raw readcount matrix. When used
#' within the Shiny app, this will be assay(SCEsetObject, "counts").
#' @param totalReads Numeric. The total number of reads in the simulated
#' dataset, to be split between all simulated cells.
#' @param cells Numeric. The number of virtual cells to simulate.
#' @param realLabels Factor. The condition labels for differential expression.
#' If only two factors present, will default to t-test. If multiple factors,
#' will default to ANOVA.
#'
#' @return A simulated counts matrix, the first row of which contains the 'true'
#' labels for each virtual cell.
#' @export
#'
generateSimulatedData <- function(totalReads, cells, originalData, realLabels){
  cells <- sample(dim(originalData)[2], size = cells, replace = TRUE)
  totalReads <- floor(totalReads / length(cells))
  originalData <- t(t(originalData) / apply(originalData, 2, sum))
  output <- matrix(nrow = dim(originalData)[1], ncol = length(cells))
  for(i in 1:length(cells)){
    output[, i] <- rmultinom(1, totalReads, originalData[, cells[i]])
  }
  return(rbind(realLabels[cells], output))
}

#' Passes the output of generateSimulatedData() to differential expression
#' tests, picking either t-tests or ANOVA for data with only two conditions or
#' multiple conditions, respectively.
#'
#' @param tempData Matrix. The output of generateSimulatedData(), where the
#' first row contains condition labels.
#'
#' @return A vector of fdr-adjusted p-values for all genes. Nonviable results
#' (such as for genes with 0 counts in a simulated dataset) are coerced to 1.
#' @export
#'
subDiffEx <- function(tempData){
  realLabels <- tempData[1, ]
  output <- tempData[-1, ]
  if(length(levels(as.factor(realLabels))) > 2){
    fdr <- subDiffEx_anova(output, realLabels)
  } else if(length(levels(as.factor(realLabels))) == 2){
    fdr <- subDiffEx_ttest(output, realLabels)
  }
  else{
    stop("Only 1 (or 0?) factor in ", levels(as.factor(realLabels)))
  }
  fdr[which(is.na(fdr))] <- 1
  return(fdr)
}

#' Returns significance data from a snapshot.
#'
#' @param originalData TODO: document
#' @param realLabels TODO: document
#' @param totalReads TODO: document
#' @param cells TODO: document
#' @param iterations TODO: document
#'
#' @return A matrix of significance information from a snapshot
#' @export
#'
iterateSimulations <- function(originalData, realLabels, totalReads, cells,
                               iterations){
  realLabels <- colData(originalData)[, realLabels]
  originalData <- counts(originalData)
  sigMatrix <- matrix(nrow = dim(originalData)[1])
  for(i in 1:iterations){
    tempData <- generateSimulatedData(totalReads, cells, originalData,
                                      realLabels = as.factor(realLabels))
    sigMatrix <- cbind(sigMatrix, subDiffEx(tempData))
  }
  return(sigMatrix[, -1])
}

#' Runs t-tests on all genes in a simulated dataset with 2 conditions, and
#' adjusts for FDR.
#'
#' @param dataset Matrix. A simulated counts matrix, sans labels.
#' @param class.labels Factor. The condition labels for the simulated cells.
#' Will be coerced into 1's and 0's.
#' @param test.type Type of test to perform. The default is t.equalvar.
#'
#' @return A vector of fdr-adjusted p-values for all genes. Nonviable results
#' (such as for genes with 0 counts in a simulated dataset) are coerced to 1.
#' @export
#'
subDiffEx_ttest <- function(dataset, class.labels, test.type = "t.equalvar") {
  class.labels <- as.numeric(as.factor(class.labels))
  class.labels <- class.labels - 1
  class.labels[class.labels > 0] <- 1
  tval <- multtest::mt.teststat(dataset, classlabel = class.labels,
                                test = test.type, nonpara = "n")
  df <- (ncol(dataset) - 2)
  pval <- 2 * (1 - pt(abs(tval), df))
  fdr <- p.adjust(pval, method = "fdr")
  return(fdr)
}

#' Runs ANOVA on all genes in a simulated dataset with more than 2 conditions,
#' and adjusts for FDR.
#'
#' @param countMatrix Matrix. A simulated counts matrix, sans labels.
#' @param condition Factor. The condition labels for the simulated cells.
#'
#' @return A vector of fdr-adjusted p-values for all genes. Nonviable results
#' (such as for genes with 0 counts in a simulated dataset) are coerced to 1.
#' @export
#'
subDiffEx_anova <- function(countMatrix, condition){
  mod <- model.matrix(~as.factor(condition))
  mod0 <- model.matrix(~1, data = condition)
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
  fstats <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1))
  p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  return(p.adjust(p, method = "fdr"))
}

#' Finds the effect sizes for all genes in the original dataset, regardless of
#' significance.
#' @param countMatrix Matrix. A simulated counts matrix, sans labels.
#' @param condition Factor. The condition labels for the simulated cells. If
#' more than 2 conditions are given, the first will be compared to all others by
#' default.
#' @return A vector of cohen's d effect sizes for each gene.
#' @export
#'
calcEffectSizes <- function(countMatrix, condition){
  groups <- levels(as.factor(unlist(condition)))
  return((apply(countMatrix[, condition == unlist(groups)[1]], 1, mean) -
            apply(countMatrix[, condition != unlist(groups)[1]], 1, mean)) /
           apply(countMatrix, 1, sd))
}

powerCalc <- function(datamatrix, sampleSizeRange=c(1000, 10000000),
                      byBatch=FALSE, batch=NULL, numSize=25) {
  if (byBatch == FALSE){
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                numSize))
    for (j in 1:dim(datamatrix)[2]) {
      probs <- as.numeric(datamatrix[, j] / sum(datamatrix[, j]))
      for (i in 1:length(probs)){
        discoveryPower <- 1 - dbinom(0, size = floor(seq.int(
          from = sampleSizeRange[1],to = sampleSizeRange[2],
          length.out = numSize)), prob = probs[i])
        outmat[i, j, ] <- discoveryPower
      }
    }
  }
  else {
    outmat <- array(NA, dim = c(dim(datamatrix)[1], dim(datamatrix)[2],
                                numSize))
    for (j in 1:nlevels(batch)) {
      probs <- datamatrix[, which(batch == levels(batch)[j])] /
        sum(datamatrix[, which(batch == levels(batch)[j])])
      for (i in 1:length(as.vector(probs))) {
        discoveryPower <- 1 - dbinom(0, size = floor(seq.int(
          from = sampleSizeRange[1], to = sampleSizeRange[2],
          length.out = numSize)), prob = as.vector(probs)[i])
        outmat[((i - 1) %% dim(datamatrix)[1]) + 1, which(batch == levels(batch)[j])[ceiling(i / dim(datamatrix)[1])], ] <- discoveryPower
      }
    }
  }
  return(outmat)
}

#' Calculate power to discover differentially expressed genes upon subsampling
#'
#' @param datamatrix The original raw readcount matrix. When used within the
#' Shiny app, this will be assay(SCEsetObject, "counts").
#' @param downmatrix The 4-dimensional array of simulated subsampled datasets.
#' Produced by Downsample()
#' @param conditions A vector of conditions, such as treatment or batch. Passed
#' as colData(SCtkExperiment)$treatment.
#' @param genelist Optional vector of genes to be searched for differential
#' expression. If provided, diffexp will not be performed on the original data.
#' @param significance Significance threshold, a float between 0 and 1
#' (exclusive).
#' @param method Which method should be used for differential expression.
#' Default is simple tpm-based t-test for speed.
#'
#' @return A matrix of recapitulation - rows are genes that were differentially
#' expressed in the original set (or passed as 'genelist'), columns are
#' simulated depths.
#' @export
#'
differentialPower <- function(datamatrix, downmatrix, conditions,
                              genelist=FALSE, significance=0.05,
                              method="tpm.t"){
  condition <- as.factor(conditions)
  if (genelist == FALSE){
    if (method == "tpm.t"){
      shrunkDown <- datamatrix[apply(datamatrix, 1, sum) > 10, ]
      genenames <- rownames(datamatrix[apply(datamatrix, 1, sum) > 10, ])
      scaledMatrix <- apply(shrunkDown, 2, function(x){
        x / sum(x)
      })
      t.result <- multiple.t.test(condition, scaledMatrix)
      genelist <- genenames[which(t.result <= significance)]
    }
    else{
      shrunkDown <- datamatrix[apply(datamatrix, 1, sum) > 10, ]
      genenames <- rownames(datamatrix[apply(datamatrix, 1, sum) > 10, ])
      countData <- newCountDataSet(shrunkDown, condition)
      countData <- estimateSizeFactors(countData)
      countData <- estimateDispersions(countData, method = "pooled",
                                       fitType = "local")
      diff.results <- nbinomTest(countData, levels(condition)[1],
                                 levels(condition)[2])
      top.results <- p.adjust(diff.results$pval, method = "fdr")
      genelist <- genenames[which(top.results <= significance)]
    }
  }
  #Create an empty matrix to keep track of how often significant genes are
  #rediscovered after downsampling
  rediscovered <- matrix(rep(0, length(genelist) * dim(downmatrix)[3]),
                         nrow = length(genelist))
  rownames(rediscovered) <- genelist
  #For each downsampling level:
  for (k in 1:dim(downmatrix)[3]) {
    #For each iteration:
    for (l in 1:dim(downmatrix)[4]) {
      if (method == "tpm.t"){
        shrunkDown <- downmatrix[apply(downmatrix[, , k, l], 1, sum) > 0, , k, l]
        genenames <- rownames(datamatrix[apply(downmatrix[, , k, l], 1, sum) > 0, ])
        scaledMatrix <- apply(shrunkDown, 2, function(x){x / sum(x)})
        t.result <- multiple.t.test(condition, scaledMatrix)
        newgenes <- genenames[which(t.result <= significance)]
        rediscovered[rownames(rediscovered) %in% newgenes, k] <- rediscovered[rownames(rediscovered) %in% newgenes, k] + 1
      }
      else{
        shrunkDown <- downmatrix[apply(downmatrix[, , k, l], 1, sum) > 0, , k, l]
        genenames <- rownames(datamatrix[apply(downmatrix[, , k, l], 1, sum) > 0, ])
        countData <- newCountDataSet(shrunkDown, condition)
        countData <- estimateSizeFactors(countData)
        countData <- estimateDispersions(countData, method = "pooled",
                                         fitType = "local")
        diff.results <- nbinomTest(countData, levels(condition)[1],
                                   levels(condition)[2])
        top.results <- p.adjust(diff.results$pval, method = "fdr")
        newgenes <- genenames[which(top.results <= significance)]
        rediscovered[rownames(rediscovered) %in% newgenes, k] <- rediscovered[rownames(rediscovered) %in% newgenes, k] + 1
      }
    }
  }
  rediscovered <- rediscovered / dim(downmatrix)[4]
  return(rediscovered)
}
