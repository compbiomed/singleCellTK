.readMatrix <- function(path) {
  
  res <- readRDS(path)
  matrix <- t(as.matrix(res))
  
  
  if (class == "Matrix") {
    
    return(matrix)
    
  } else if (class == "DelayedArray") {
    
    res <- DelayedArray::DelayedArray(res)
    
    return(matrix)
  }
    
    
  
}


.readMetrics <- function(path) {
  
  res <- fread(path)
  res <- as.dataframe(res)
  return(res)
  
  if(ncol(res)==1)
  {
    stop("There are no Cell or Gene Metrics!")
  }
  
}

.readEmptyDrops<-function(path){
  EmptyDrops<-read.csv(path)
}


.constructSCEFromOptimusOutputs <- function(matrixLocation,
                                               
                                               CellMetricsLocation,
                                               
                                               EmptyDropsLocation,
                                               
                                               GeneMetricsLocation) {
  
  
  c_me <- .readMetrics(file.path(CellMetricsLocation))
  
  cb <- c_me[,1]
  
  c_me <- c_me[,-1]
  
  emptydrops <- .readEmptyDrops(file.path(EmptyDropsLocation))
  
  c_me_full <- cbind(c_me, emptydrops)
  
  g_me <- .readMetrics(file.path(GeneMetricsLocation))
  
  g_me <- g_me[,g_me %in% rownames(ma)]
  
  gi <- g_me[,1]
  
  g_me <- g_me[,-1]
  
  ma <- .readMatrix(file.path(matrixLocation))
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    
    assays = list(counts = ma))
  
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(g_me,row.names = gi)
  
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(c_me_full,row.names = cb)
  
  
  
  return(sce)
  
}


.checkArgsImportOptimus <- function(OptimusDirs) {
  
  if (!dir.exists(OptimusDirs)) {
      stop("OptimusDirs is NULL!")
  }
}

.importOptimus <- function(
  
  OptimusDirs,
  
  matrixLocation,
  
  CellMetricsLocation,
  
  EmptyDropsLocation,
  
  GeneMetricsLocation) {
  
  
  .checkArgsImportOptimus(OptimusDirs)
    
  sce <- .constructSCEFromOptimusOutputs(OptimusDirs,
                                              
                                             matrixLocation = matrixLocation,
                                               
                                             CellMetricsLocation = CellMetricsLocation,
                                         
                                             EmptyDropsLocation = EmptyDropsLocation,
                                               
                                             GeneMetricsLocation = GeneMetricsLocation)
    
  
   return(sce)
  
}


