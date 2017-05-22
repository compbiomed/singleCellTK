#' thresholdGenes 
#' 
#' @param SCEdata SCESet object
#' @export
thresholdGenes <- function(SCEdata){
  # data preparation 
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, pData(SCEdata), fdata)
  
  SCE_new <- SCE_new[which(freq(SCE_new)>0),]
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  return(thres)
}

#' MAST Violin
#' 
#' Run MAST analysis on a SCESet object.
#' 
#' @param SCEdata SCESet object
#' @param fcHurdleSig The filtered result from hurdle model 
#' @param samplesize The number of most significant genes
#' @param threshP Plot threshold values from adaptive thresholding. Default is
#' FALSE
#' @param variable Select the condition of interest
#' @export
MASTviolin <- function(SCEdata, fcHurdleSig, samplesize = 49, threshP=FALSE,
                       variable){
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, pData(SCEdata), fdata)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new)>0),]
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                     min_per_bin = 30)
  assays(SCE_new) <- list(thresh=thres$counts_threshold, tpm=assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:49]
  flat_dat <- as(SCE_new[entrez_to_plot,], 'data.table')
  if(threshP){
    ggbase <- ggplot(flat_dat,
                     aes_string(x=variable, y="thresh", color=variable)) + 
      geom_jitter() + facet_wrap(~primerid, scale='free_y', ncol=7)
  }
  else{
    ggbase <- ggplot(flat_dat,
                     aes_string(x=variable, y="tpm", color=variable)) + 
      geom_jitter() + facet_wrap(~primerid, scale='free_y', ncol=7)
  }
  violinplot <- ggbase + geom_violin() + ggtitle("Violin Plot")
  return(violinplot)
}

#' MAST linear model
#' 
#' Run MAST analysis on a SCESet object.
#' 
#' @param SCEdata SCESet object
#' @param fcHurdleSig The filtered result from hurdle model 
#' @param samplesize The number of most significant genes
#' @export
MASTregression <- function(SCEdata,fcHurdleSig, 
                           samplesize = 49){
  count <- counts(SCEdata)
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  SCE_new <- MAST::FromMatrix(expres, pData(SCEdata), fdata)
  SCE_new_lm <- MAST::FromMatrix(count, pData(SCEdata), fdata)
  cdr2 <-colSums(assay(SCE_new_lm)>0)
  colData(SCE_new)$cngeneson <- scale(cdr2)
  SCE_new <- SCE_new[which(MAST::freq(SCE_new)>0),]
  
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20, min_per_bin = 30)
  assays(SCE_new) <- list(thresh=thres$counts_threshold, tpm=assay(SCE_new))
  entrez_to_plot <- fcHurdleSig$Gene[1:49]
  flat_dat <- as(SCE_new[entrez_to_plot,], 'data.table')
  ggbase <- ggplot(flat_dat, aes(x=level1class, y=tpm,color=level1class)) + geom_jitter()+facet_wrap(~primerid, scale='free_y',ncol=7)
  flat_dat[,lmPred:=lm(thresh~cngeneson + level1class)$fitted, key=primerid]
  regressionplot <- ggbase +aes(x=cngeneson) + geom_line(aes(y=lmPred), lty=1) + xlab('Standardized Cellular Detection Rate')
  return(regressionplot)
}
