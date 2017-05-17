## What the hell is the $ operator is invalid for atomic vectors and subscript out of bounds


#' MAST 
#' @param data SCE dataset for running the MAST 
#' @param FCTHRESHOLD the threshold for filetering the sample 
#' @param freq_expressed 
#' @example TODO
#' @export

MAST <- function(SCEdata, 
                 lbound = 0.1, 
                 FCTHRESHOLD=log2(1.5), freq_expressed= 0.2 ){
  library(data.table)
  library(MAST)
  # data preparation 
  nGeneOn <- colSums(counts(SCEdata))
  pdata <- pData(SCEdata)
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  count <- counts(SCEdata)
  SCE_new <- FromMatrix(count, pdata, fdata)
  
  # filter
  SCE_new_sample <- SCE_new[sample(which(freq(SCE_new)>lbound)),]
  
  
  # Defferential expression using a hurdle model using the filtered dataset
  hurdle1 <- zlm.SingleCellAssay(~level1class,SCE_new_sample)
  summaryh1 <- MAST::summary(hurdle1,doLRT='level1classmicroglia')
  
  summaryDT <- summaryh1[['datatable']]
  
  # 
  fcHurdle <- merge(summaryDT[summaryDT$contrast=='level1classmicroglia' & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],
                    summaryDT[summaryDT$contrast=='level1classmicroglia' & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')])
  #fcHurdle <- data.table(fcHurdle)
  
  # Use p-value correction method, here we use fdr
  fcHurdle$fdr <- p.adjust(fcHurdle[,"Pr(>Chisq)"],'fdr')
  resultdata <- fcHurdle
  
  # Filter the data again by the adjusted pvalue and coef
  fcHurdleSig <-  merge(fcHurdle[fcHurdle$fdr<.05 & abs(fcHurdle$coef)>FCTHRESHOLD,], 
                        as.data.table(mcols(SCE_new_sample)), by='primerid')
  setorder(fcHurdleSig, fdr)
  
  return(fcHurdleSig)
  
}