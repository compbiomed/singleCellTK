#' Downsample Data
#'
#' @param datamatrix TODO:document
#' @param newcounts TODO:document
#' @param byBatch TODO:document
#' @param batch TODO:document
#' @param iterations TODO:document
#'
#' @return Downsampled matrix
#' @export Downsample
Downsample <- function(datamatrix, newcounts = c(4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144),
                              byBatch = FALSE, batch = NULL, iterations = 10) {
  if (byBatch == FALSE) {
    outmat <- array(, dim = c(dim(datamatrix)[1], dim(datamatrix)[2], length(newcounts), iterations))
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
    outmat <- array(, dim = c(dim(datamatrix)[1], dim(datamatrix)[2], length(newcounts), iterations))
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

generateSimulatedData <- function(totalReads, cells, originalData, realLabels){
  cells <- sample(dim(originalData)[2], size=cells, replace=TRUE)
  totalReads <- floor(totalReads/length(cells))
  originalData <- t(t(originalData)/apply(originalData,2,sum))
  output <- matrix(nrow=dim(originalData)[1],ncol=length(cells))
  for(i in 1:length(cells)){
    output[,i] <- rmultinom(1,totalReads,originalData[,cells[i]])
  }
  if(length(levels(as.factor(realLabels))) > 2){
    fdr <- subDiffEx_anova(output,realLabels[cells])
  }
  else{
    fdr <- subDiffEx_ttest(output,realLabels[cells])
  }
  fdr[which(is.na(fdr))] <- 1
  return(fdr)
}

iterateSimulations <- function(originalData, realLabels, totalReads, cells, iterations){
  sigMatrix <- matrix(nrow=dim(originalData)[1])
  for(i in 2:iterations){
    sigMatrix <- cbind(sigMatrix, generateSimulatedData(totalReads, cells, originalData, realLabels=as.factor(realLabels)))
  }
  return(sigMatrix)
}

require(multtest)
subDiffEx_ttest <- function(dataset, class.labels, test.type = "t.equalvar") {
  tval <- multtest::mt.teststat(dataset, classlabel = class.labels, test = test.type, nonpara = "n")
  df <- (ncol(dataset) - 2)
  pval <- 2 * (1 - pt(abs(tval), df))
  fdr <- p.adjust(pval, method = "fdr")
  return(fdr)
}

subDiffEx_anova <- function(countMatrix, condition){
  mod <- model.matrix(~as.factor(condition))
  mod0 <- model.matrix(~1,data=condition)
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
  rss1 <- resid^2 %*% rep(1, n)
  rss0 <- resid0^2 %*% rep(1, n)
  fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
  p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  return(p.adjust(p, method = "fdr"))
}


calcEffectSizes <- function(countMatrix, condition){
  groups <- levels(as.factor(condition))
  return((apply(countMatrix[,condition==groups[1]],1,mean)-apply(countMatrix[,condition!=groups[1]],1,mean))/apply(countMatrix,1,sd))
}

powerCalc <- function(datamatrix, sampleSizeRange=c(1000, 10000000), byBatch=FALSE, batch=NULL, numSize=25) {
  if (byBatch == FALSE){
    outmat <- array(, dim = c(dim(datamatrix)[1], dim(datamatrix)[2], numSize))
    for (j in 1:dim(datamatrix)[2]) {
      probs <- as.numeric(datamatrix[, j] / sum(datamatrix[, j]))
      for (i in 1:length(probs)){
        discoveryPower <- 1 - dbinom(0, size = floor(seq.int(from = sampleSizeRange[1], to = sampleSizeRange[2], length.out = numSize)), prob = probs[i])
        outmat[i, j, ] <- discoveryPower
      }
    }
  }
  else {
    outmat <- array(, dim = c(dim(datamatrix)[1], dim(datamatrix)[2], numSize))
    for (j in 1:nlevels(batch)) {
      probs <- datamatrix[, which(batch == levels(batch)[j])] / sum(datamatrix[, which(batch == levels(batch)[j])])
      for (i in 1:length(as.vector(probs))) {
        discoveryPower <- 1 - dbinom(0, size = floor(seq.int(from = sampleSizeRange[1], to = sampleSizeRange[2], length.out = numSize)), prob = as.vector(probs)[i])
        outmat[((i - 1) %% dim(datamatrix)[1]) + 1, which(batch == levels(batch)[j])[ceiling(i / dim(datamatrix)[1])], ] <- discoveryPower
      }
    }
  }
  return(outmat)
}


#' Calculate power to discover differentially expressed genes upon subsampling
#'
#' @param datamatrix The original raw readcount matrix. When used within the Shiny app, this will be assay(SCEsetObject, "counts").
#' @param downmatrix The 4-dimensional array of simulated subsampled datasets. Produced by Downsample()
#' @param conditions A vector of conditions, such as treatment or batch. Passed as colData(SingleCellExperiment)$treatment.
#' @param genelist Optional vector of genes to be searched for differential expression. If provided, diffexp will not be performed on the oringinal data.
#' @param significance Significance threshold, a float between 0 and 1 (exclusive).
#' @param method Which method should be used for differential expression. Default is simple tpm-based t-test for speed.
#'
#' @return A matrix of recapitulation - rows are genes that were differentially expressed in the original set (or passed as 'genelist'), columns are simulated depths.
#' @export differentialPower
differentialPower <- function(datamatrix, downmatrix, conditions, genelist=FALSE, significance=0.05, method="tpm.t") {
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
      countData <- estimateDispersions(countData, method = "pooled", fitType = "local")
      diff.results <- nbinomTest(countData, levels(condition)[1], levels(condition)[2])
      top.results <- p.adjust(diff.results$pval, method = "fdr")
      genelist <- genenames[which(top.results <= significance)]
    }
  }
  #Create an empty matrix to keep track of how often significant genes are rediscovered after downsampling
  rediscovered <- matrix(rep(0, length(genelist) * dim(downmatrix)[3]), nrow = length(genelist))
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
        countData <- estimateDispersions(countData, method = "pooled", fitType = "local")
        diff.results <- nbinomTest(countData, levels(condition)[1], levels(condition)[2])
        top.results <- p.adjust(diff.results$pval, method = "fdr")
        newgenes <- genenames[which(top.results <= significance)]
        rediscovered[rownames(rediscovered) %in% newgenes, k] <- rediscovered[rownames(rediscovered) %in% newgenes, k] + 1
      }
    }
  }
  rediscovered <- rediscovered / dim(downmatrix)[4]
  return(rediscovered)
}
