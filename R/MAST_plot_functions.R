#' thresholdGenes 
#' 
#' @param SCEdata SCESet object
#' @export
thresholdGenes <- function(SCEdata){
  # data preparation 
  expres <- exprs(SCEdata)
  fdata = data.frame(Gene = rownames(expres))
  rownames(fdata) = fdata$Gene
  count <- counts(SCEdata)
  SCE_new <- MAST::FromMatrix(expres, pData(SCEdata), fdata)
  
  SCE_new <- SCE_new[which(freq(SCE_new)>0),]
  thres <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20, min_per_bin = 30)
  return(thres)
}
