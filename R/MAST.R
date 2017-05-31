#' MAST 
#' 
#' Run MAST analysis on a SCESet object.
#' 
#' @param SCEdata SCESet object
#' @param condition select varible (from the pdata) that is used for modle
#' @param interest.level If the condition of interest has more than two factors,
#' indicate which level should be used to compare to all other samples.
#' @param freq_expressed Filter genes that are expressed in at least this
#' fraction of cells. The default is expression in 0.1 of samples.
#' @param fc_threshold Minimum fold change for differentially expressed gene.
#' @param p.value p values for selecting the hurdle result, default is 0.05
#' @param usethresh Use adaptive thresholding to filter genes. The default is
#' FALSE.
#' @export
MAST <- function(SCEdata, condition = NULL, interest.level = NULL,
                 freq_expressed = 0.1, fc_threshold=log2(1.5), p.value = 0.05,
                 usethresh=FALSE){
  
  if(is.null(condition)){
    stop("specify the condition of interest")
  }
  
  if(length(unique(pData(SCEdata)[,condition])) == 1) {
    stop("only one level is in the condition")
  }
  
  if(is.null(interest.level) & length(unique(pData(SCEdata)[, condition]))>2 ){
    stop("You must specify a level of interest when more than 2 levels are in",
         " the condition")
  }

  # Create MAST SingleCellAssay 
  pdata <- pData(SCEdata)
  expres <- exprs(SCEdata)
  fdata <- fData(SCEdata)
  SCE_new <- MAST::FromMatrix(expres, pdata, fdata)
  
  #Caculate CDR for zlm model
  colData(SCE_new)$cngeneson <- scale(colSums(assay(SCE_new) > 0))
  
  if(usethresh){
    SCE_new <- SCE_new[which(freq(SCE_new)>0),]
    thresh <- thresholdSCRNACountMatrix(assay(SCE_new), nbins = 20,
                                       min_per_bin = 30)
    assays(SCE_new) <- list(thresh=thresh$counts_threshold, tpm=assay(SCE_new))
  }
  
  # filter based on frequency of expression across samples
  SCE_new_sample <- SCE_new[which(MAST::freq(SCE_new) > freq_expressed),]
  
  # if the condition of interest is numeric, to change it to a factor
  if(is.numeric(colData(SCE_new_sample)[,condition])){  
    colData(SCE_new_sample)[,condition] <- as.factor(colData(SCE_new_sample)[,condition])
  }
  
  # >2 levels in the condition 
  if(!is.null(interest.level) & length(unique(pData(SCEdata)[, condition])) > 2){
    levels(colData(SCE_new_sample)[,condition]) <- c(levels(colData(SCE_new_sample)[,condition]), paste0("no_", interest.level))
    colData(SCE_new_sample)[,condition][colData(SCE_new_sample)[,condition] != interest.level] <- paste0("no_", interest.level)
    colData(SCE_new_sample)[,condition] <- droplevels(as.factor(colData(SCE_new_sample)[,condition]))
    
    hurdle1 <- MAST::zlm(as.formula(paste0('~', condition, "+cngeneson")),SCE_new_sample)
    summaryh1 <- MAST::summary(hurdle1,doLRT=paste0(condition, "no_", interest.level))  
    
    summaryDT <- summaryh1[['datatable']]
    
    fcHurdle <- merge(summaryDT[summaryDT$contrast==paste0(condition, "no_", interest.level) & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],   # hurdel p value 
                      summaryDT[summaryDT$contrast==paste0(condition, "no_", interest.level) & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')]  # logFC coefficients
    )
  } else {
    colData(SCE_new_sample)[,condition] <- droplevels(as.factor(colData(SCE_new_sample)[,condition]))
    level.cond  <- levels(colData(SCE_new_sample)[,condition])
    
    hurdle1 <- MAST::zlm(as.formula(paste0('~', condition, "+cngeneson")),SCE_new_sample)
    summaryh1 <- MAST::summary(hurdle1,doLRT=paste0(condition, level.cond[2]))  
    
    summaryDT <- summaryh1[['datatable']]
    
    fcHurdle <- merge(summaryDT[summaryDT$contrast==paste0(condition, level.cond[2]) & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],   # hurdel p value 
                      summaryDT[summaryDT$contrast==paste0(condition, level.cond[2]) & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')]  # logFC coefficients
    )
  }
  
  # Use p-value correction method, here we use fdr
  fcHurdle$fdr <- p.adjust(fcHurdle$'Pr(>Chisq)', 'fdr')
  
  # Filter the data again by the adjusted pvalue and coef
  fcHurdleSig <-  fcHurdle[fcHurdle$fdr<p.value & abs(fcHurdle$coef) > fc_threshold & !is.nan(fcHurdle$coef),]
  colnames(fcHurdleSig)[1] <- "Gene"
                      
  data.table::setorder(fcHurdleSig, fdr)
  
  return(fcHurdleSig)
}
