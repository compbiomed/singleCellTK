#' MAST 
#' 
#' Run MAST analysis on a SCESet object.
#' 
#' @param SCEdata SCESet object
#' @param FCTHRESHOLD the threshold for filetering the sample 
#' @param freq_expressed Filter genes that are expressed in at least this fraction of 
#' cells
#' @param p.value p-value cutoff
#' @export
MAST <- function(SCEdata, 
                 FCTHRESHOLD=log2(1.5), freq_expressed= 0.2,
                 p.value=0.05){
  # data preparation 
  nGeneOn <- colSums(counts(SCEdata))
  pdata <- pData(SCEdata)
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  count <- counts(SCEdata)
  SCE_new <- MAST::FromMatrix(count, pdata, fdata)
  
  # filter
  SCE_new_sample <- SCE_new[sample(which(MAST::freq(SCE_new)>freq_expressed)),]
  
  
  # Defferential expression using a hurdle model using the filtered dataset
  hurdle1 <- MAST::zlm(~level1class,SCE_new_sample)
  summaryh1 <- MAST::summary(hurdle1,doLRT='level1classoligodendrocytes')
  
  summaryDT <- summaryh1[['datatable']]
  
  # 
  fcHurdle <- merge(summaryDT[summaryDT$contrast=='level1classoligodendrocytes' & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],
                    summaryDT[summaryDT$contrast=='level1classoligodendrocytes' & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')])
  
  # Use p-value correction method, here we use fdr
  fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`,'fdr')
  resultdata <- fcHurdle
  
  # Filter the data again by the adjusted pvalue and coef
  fcHurdleSig <-  fcHurdle[fcHurdle$fdr<p.value & abs(fcHurdle$coef)>FCTHRESHOLD,]
  colnames(fcHurdleSig)[1] <- "Gene"
  data.table::setorder(fcHurdleSig, fdr)
  
  return(fcHurdleSig)
}
