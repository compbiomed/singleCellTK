#' Downsample Data
#'
#' @param datamatrix 
#' @param newcounts 
#' @param byBatch 
#' @param batch 
#' @param iterations 
#'
#' @return Downsampled matrix
#' @export Downsample
Downsample <- function(datamatrix, newcounts = c(4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144),
                              byBatch = FALSE, batch = NULL, iterations = 10) {
  if (byBatch == FALSE) {
    outmat <- array(, dim=c(dim(datamatrix)[1], dim(datamatrix)[2], length(newcounts), iterations))
    for (j in 1:dim(datamatrix)[2]) {
      probs <- datamatrix[, j] / sum(datamatrix[, j])
      for (k in 1:length(newcounts)) {
        samps <- rmultinom(iterations, newcounts[k], probs)
        for (l in 1:iterations) {
          outmat[,j,k,l] <- samps[,l]
        }
      }
    }
  }
  else {
    outmat <- array(, dim=c(dim(datamatrix)[1], dim(datamatrix)[2], length(newcounts), iterations))
    for (j in 1:nlevels(batch)) {
      probs <- datamatrix[,which(batch == levels(batch)[j])] / sum(datamatrix[,which(batch == levels(batch)[j])])
      for (k in 1:length(newcounts)) {
        samps <- rmultinom(iterations, newcounts[k], as.vector(probs))
        for (l in 1:iterations) {
          outmat[,which(batch == levels(batch)[j]),k,l] <- as.matrix(samps[,l], nrow = dim(datamatrix)[1])
        }
      }
    }
  }
  return(outmat)
}

powerCalc <- function(datamatrix, sampleSizeRange=c(1000,10000000), byBatch=FALSE, batch=NULL, numSize=25) {
  if (byBatch == FALSE){
    outmat <- array(, dim=c(dim(datamatrix)[1], dim(datamatrix)[2], numSize))
    for (j in 1:dim(datamatrix)[2]) {
      probs <- as.numeric(datamatrix[, j] / sum(datamatrix[, j]))
      for(i in 1:length(probs)){
        discoveryPower <- 1 - dbinom(0, size=floor(seq.int(from=sampleSizeRange[1], to=sampleSizeRange[2], length.out=numSize)), prob=probs[i])
        outmat[i,j,] <- discoveryPower
      }
    }    
  }
  else {
    outmat <- array(, dim=c(dim(datamatrix)[1], dim(datamatrix)[2], numSize))
    for (j in 1:nlevels(batch)) {
      probs <- datamatrix[,which(batch == levels(batch)[j])] / sum(datamatrix[,which(batch == levels(batch)[j])])
      for(i in 1:length(as.vector(probs))) {
        discoveryPower <- 1-dbinom(0, size=floor(seq.int(from=sampleSizeRange[1], to=sampleSizeRange[2], length.out=numSize)), prob=as.vector(probs)[i])
        outmat[((i-1)%%dim(datamatrix)[1])+1,which(batch == levels(batch)[j])[ceiling(i/dim(datamatrix)[1])],] <- discoveryPower
      }
    }
  }
  return(outmat)
}

require(multtest)
multiple.t.test <- function(class.labels, dataset, test.type="t.equalvar") {
  tval <- multtest::mt.teststat(dataset, classlabel=class.labels, test=test.type, nonpara="n")
  df <- (ncol(dataset) - 2)
  pval <- 2*(1-pt(abs(tval),df))
  fdr <- p.adjust(pval, method="fdr")
  return(fdr)
}
require(DESeq)

#' Calculate power to discover differentially expressed genes upon subsampling
#'
#' @param datamatrix The original raw readcount matrix. When used within the Shiny app, this will be counts(SCEsetObject).
#' @param downmatrix The 4-dimensional array of simulated subsampled datasets. Produced by Downsample()
#' @param conditions A vector of conditions, such as treatment or batch. Passed as phenoData(SCEsetObject)$treatment.
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
      shrunkDown <- datamatrix[apply(datamatrix,1,sum) > 10, ]
      genenames <- rownames(datamatrix[apply(datamatrix,1,sum) > 10, ])
      scaledMatrix <- apply(shrunkDown,2,function(x){x/sum(x)})
      t.result <- multiple.t.test(condition, scaledMatrix)
      genelist <- genenames[which(t.result <= significance)]
    }
    else{
      shrunkDown <- datamatrix[apply(datamatrix,1,sum) > 10, ]
      genenames <- rownames(datamatrix[apply(datamatrix,1,sum) > 10, ])
      countData <- newCountDataSet( shrunkDown, condition )
      countData <- estimateSizeFactors( countData )
      countData <- estimateDispersions( countData, method="pooled", fitType="local" )
      diff.results <- nbinomTest( countData, levels(condition)[1], levels(condition)[2] )
      top.results <- p.adjust( diff.results$pval, method="fdr" )
      genelist <- genenames[which(top.results <= significance)]
    }
  }
  #Create an empty matrix to keep track of how often significant genes are rediscovered after downsampling
  rediscovered <- matrix(rep(0, length(genelist)*dim(downmatrix)[3]), nrow=length(genelist))
  rownames(rediscovered) <- genelist
  #For each downsampling level:
  for (k in 1:dim(downmatrix)[3]) {
    #For each iteration:
    for (l in 1:dim(downmatrix)[4]) {
      if (method == "tpm.t"){
        shrunkDown <- downmatrix[apply(downmatrix[,,k,l],1,sum) > 0, ,k, l]
        genenames <- rownames(datamatrix[apply(downmatrix[,,k,l],1,sum) > 0, ])
        scaledMatrix <- apply(shrunkDown,2,function(x){x/sum(x)})
        t.result <- multiple.t.test(condition, scaledMatrix)
        newgenes <- genenames[which(t.result <= significance)]
        rediscovered[rownames(rediscovered) %in% newgenes, k] <- rediscovered[rownames(rediscovered) %in% newgenes, k] + 1
      }
      else{
        shrunkDown <- downmatrix[apply(downmatrix[,,k,l],1,sum) > 0, ,k, l]
        genenames <- rownames(datamatrix[apply(downmatrix[,,k,l],1,sum) > 0, ])
        countData <- newCountDataSet( shrunkDown, condition )
        countData <- estimateSizeFactors( countData )
        countData <- estimateDispersions( countData, method="pooled", fitType="local" )
        diff.results <- nbinomTest( countData, levels(condition)[1], levels(condition)[2])
        top.results <- p.adjust( diff.results$pval, method="fdr" )
        newgenes <- genenames[which(top.results <= significance)]
        rediscovered[rownames(rediscovered) %in% newgenes, k] <- rediscovered[rownames(rediscovered) %in% newgenes, k] + 1
      }
    }
  }
  rediscovered <- rediscovered/dim(downmatrix)[4]
  return(rediscovered)
}