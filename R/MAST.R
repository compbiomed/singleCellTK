#' MAST 
#' 
#' #' Run MAST analysis on a SCESet object.
#' 
#' @param SCEdata SCESet object
#' @param FCTHRESHOLD the threshold for filetering the sample 
#' @param freq_expressed Filter genes that are expressed in at least this fraction of 
#' cells
#' @param condition select varible (from the pdata) that is used for modle
#' @param interest.level interested level in the condition 
#' @param p.value  p values for selecting the hurdle reslut, default is 0.05
#' @export

MAST <- function(SCEdata, 
                 freq_expressed = 0.1, condition = NULL ,  interest.level = NULL, 
                 FCTHRESHOLD=log2(1.5) , p.value = 0.05 ){
  
  if(is.null(condition)){
    stop("specify the condition of interst")
  }
  
  if(length(unique(pData(SCEdata)[,condition])) == 1) {
    stop("only one level is in the condition")
  }
  
  if(is.null(interest.level) & length(unique(pData(SCEdata)[, condition]))>2 ){
    stop("have to specify a level of interest when more than 2 levels are in the condiciton")
  }
  
  # data preparation 
  nGeneOn <- colSums(counts(SCEdata))
  pdata <- pData(SCEdata)
  expres <- exprs(SCEdata)
  fdata <- fData(SCEdata)
  count <- counts(SCEdata)
  
  SCE_new <- MAST::FromMatrix(expres, pdata, fdata)
  
  # filter
  SCE_new_sample <- SCE_new[which(MAST::freq(SCE_new)>freq_expressed),]
  
  # Defferential expression using a hurdle model using the filtered dataset
  if(is.numeric(colData(SCE_new_sample)[,condition])) {  # if the condition is numeric, to change it to factor
    colData(SCE_new_sample)[,condition] <- as.factor(colData(SCE_new_sample)[,condition])
  }
  ##   >2 levels in the condition 
  if(!is.null(interest.level)){
    levels(MAST::colData(SCE_new_sample)[,condition]) <- c(levels(colData(SCE_new_sample)[,condition]), paste0("no_", interest.level))
    colData(SCE_new_sample)[,condition][colData(SCE_new_sample)[,condition] != interest.level] <- paste0("no_", interest.level)
    colData(SCE_new_sample)[,condition] <- droplevels(as.factor(colData(SCE_new_sample)[,condition]))
    
    #hurdle1 <- MAST::zlm.SingleCellAssay(as.formula(paste0('~', condition)),SCE_new_sample)
    hurdle1 <- MAST::zlm(as.formula(paste0('~', condition)),SCE_new_sample)
    summaryh1 <- MAST::summary(hurdle1,doLRT=paste0(condition, "no_", interest.level))  
    
    summaryDT <- summaryh1[['datatable']]
    
    fcHurdle <- merge(summaryDT[summaryDT$contrast==paste0(condition, "no_", interest.level) & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],   # hurdel p value 
                      summaryDT[summaryDT$contrast==paste0(condition, "no_", interest.level) & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')]  # logFC coefficients
    )  
  }else if(is.null(interest.level)){
    colData(SCE_new_sample)[,condition] <- droplevels(as.factor(colData(SCE_new_sample)[,condition]))
    level.cond  <- levels(colData(SCE_new_sample)[,condition])
    
    hurdle1 <- MAST::zlm(as.formula(paste0('~', condition)),SCE_new_sample)
    summaryh1 <- MAST::summary(hurdle1,doLRT=paste0(condition, level.cond[2]))  
    
    summaryDT <- summaryh1[['datatable']]
    
    fcHurdle <- merge(summaryDT[summaryDT$contrast==paste0(condition, level.cond[2]) & summaryDT$component=='H', c('primerid','Pr(>Chisq)')],   # hurdel p value 
                      summaryDT[summaryDT$contrast==paste0(condition, level.cond[2]) & summaryDT$component=='logFC', c('primerid','coef','ci.hi','ci.lo')]  # logFC coefficients
    )
  }
  
  
  
  # Use p-value correction method, here we use fdr
  fcHurdle$fdr <- p.adjust(fcHurdle$'Pr(>Chisq)','fdr')
  resultdata <- fcHurdle
  
  # Filter the data again by the adjusted pvalue and coef
  fcHurdleSig <-  fcHurdle[fcHurdle$fdr<p.value & abs(fcHurdle$coef)>FCTHRESHOLD & !is.nan(fcHurdle$coef),]
  colnames(fcHurdleSig)[1] <- "Gene"
                      
  data.table::setorder(fcHurdleSig, fdr)
  
  return(fcHurdleSig)
}
