########################################
### Helper Functions ##################
#######################################

#' .updateAssaySCEFromScanpy
#' Update/Modify/Add an assay in the provided SingleCellExperiment object from
#' an AnnData object
#' @param inSCE Input SingleCellExperiment object
#' @param scanpyObject Input annData object
#' @param assaySlotSCE Selected assay to update in the input
#' SingleCellExperiment object
#' @param scanpyAssaySlot Selected assay from annData object. Default
#' \code{"X"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  data from annData object appended to the \link{assay} slot.
#' @importFrom SummarizedExperiment assay<-
#' @noRd
.updateAssaySCEFromScanpy <- function(inSCE,
                                      scanpyObject,
                                      assaySlotSCE,
                                      scanpyAssaySlot = "X") {
  assay(inSCE, assaySlotSCE) <- NULL
  temp.matrix <- t(scanpyObject[[scanpyAssaySlot]])
  colnames(temp.matrix) <- colnames(inSCE)
  rownames(temp.matrix) <- rownames(inSCE)
  assay(inSCE, assaySlotSCE) <- temp.matrix
  
  return(inSCE)
}


#####################################################
#### Normalization and Scaling function ############
#####################################################

#' runScanpyNormalizeData
#' Wrapper for NormalizeData() function from scanpy library
#' Normalizes the sce object according to the input parameters
#' @param inSCE (sce) object to normalize
#' @param useAssay Assay containing raw counts to use for normalization.
#' @param targetSum If NULL, after normalization, each observation (cell) has a 
#' total count equal to the median of total counts for observations (cells) 
#' before normalization. Default \code{1e4}
#' @param maxFraction Include cells that have more counts than max_fraction of 
#' the original total counts in at least one cell. Default \code{0.05}
#' @param normAssayName Name of new assay containing normalized data. Default
#' \code{scanpyNormData}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' rownames(sce) <- rowData(sce)$feature_name
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyNormalizeData <- function(inSCE,
                                   useAssay,
                                   targetSum = 1e4,
                                   maxFraction = 0.05,
                                   normAssayName = "scanpyNormData") {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                             X_name = useAssay,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = FALSE,
                                             metadata = FALSE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = FALSE,
                                             verbose = NULL)
  normValue <- sc$pp$normalize_total(scanpyObject, 
                                     target_sum = targetSum, 
                                     max_fraction =  maxFraction,
                                     inplace = FALSE)
  
  scanpyObject$obs$n_counts <- normValue$norm_factor
  scanpyObject$X <- normValue$X
  
  # log transform the data.
  scanpyObject <- sc$pp$log1p(scanpyObject, copy = TRUE)
  
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, normAssayName)
  metadata(inSCE)$scanpy$normValues <- normValue$X
  colData(inSCE)$n_counts <- 
    as.factor(unlist(scanpyObject$obs['n_counts']))
  
  inSCE <- expSetDataTag(inSCE = inSCE,
                         assayType = "normalized",
                         assays = normAssayName)
  return(inSCE)
}

#' runScanpyScaleData
#' Scales the input sce object according to the input parameters
#' @param inSCE (sce) object to scale
#' @param useAssay Assay containing normalized counts to scale.
#' @param scaledAssayName Name of new assay containing scaled data. Default
#' \code{scanpyScaledData}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' }
#' @return Scaled \code{SingleCellExperiment} object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyScaleData <- function(inSCE,
                               useAssay = "scanpyNormData",
                               scaledAssayName = "scanpyScaledData") {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                             X_name = useAssay,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = FALSE,
                                             metadata = TRUE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = FALSE,
                                             verbose = NULL)
  sc$pp$scale(scanpyObject, max_value=10)
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, scaledAssayName, "X")
  rowData(inSCE)$"Scanpy_mean" <- scanpyObject$var$mean 
  rowData(inSCE)$"Scanpy_std" <- scanpyObject$var$std 
  
  inSCE <- expSetDataTag(inSCE = inSCE,
                         assayType = "scaled",
                         assays = scaledAssayName)
  
  return(inSCE)
}

######################################
### High Variable Genes ###############
######################################

#' runScanpyFindHVG
#' Find highly variable genes and store in the input sce object
#' @param inSCE (sce) object to compute highly variable genes from and to store
#' back to it
#' @param useAssay Specify the name of the assay to use for computation
#'  of variable genes. It is recommended to use log normalized data, except when 
#'  flavor='seurat_v3', in which counts data is expected.
#' @param method selected method to use for computation of highly variable
#' genes. One of \code{'seurat'}, \code{'cell_ranger'}, or \code{'seurat_v3'}.
#' Default \code{"seurat"}.
#' @param altExpName Character. Name of the alternative experiment object to
#' add if \code{returnAsAltExp = TRUE}. Default \code{featureSubset}.
#' @param altExp Logical value indicating if the input object is an
#' altExperiment. Default \code{FALSE}.
#' @param hvgNumber numeric value of how many genes to select as highly
#' variable. Default \code{2000}
#' @param minMean If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'. Default \code{0.0125}
#' @param maxMean If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'. Default \code{3}
#' @param minDisp If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'. Default \code{0.5}
#' @param maxDisp If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'. Default \code{Inf}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' g <- getTopHVG(sce, method = "seurat", hvgNumber = 500)
#' }
#' @return Updated \code{SingleCellExperiment} object with highly variable genes
#' computation stored
#' \code{\link{getTopHVG}}, \code{\link{plotTopHVG}}
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors metadata<-
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyFindHVG <- function(inSCE,
                             useAssay = "scanpyNormData",
                             method = c("seurat", "cell_ranger", "seurat_v3"),
                             altExpName = "featureSubset",
                             altExp = FALSE,
                             hvgNumber = 2000,
                             minMean = 0.0125,
                             maxMean = 3,
                             minDisp = 0.5,
                             maxDisp = Inf) {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  method <- match.arg(method)
  fullSCE_cellranger <- NULL
  
  if (!altExp) {
    if (method == "cell_ranger"){
      # for cell_ranger method only
      fullSCE_cellranger <- inSCE
      inSCE <- inSCE[which(rowSums(assay(inSCE, useAssay)) > 0), ]
      #############################
    }
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                               X_name = useAssay,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               reducedDims = FALSE,
                                               metadata = TRUE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
  }
  else{
    if (method == "cell_ranger"){
      # for cell_ranger method only
      tempSCE <- inSCE
      inSCE <- altExp(inSCE, altExpName)[which(rowSums(assay(altExp(inSCE, altExpName), useAssay)) > 0), ]
      #############################
    }
    scanpyObject <- zellkonverter::SCE2AnnData(sce = altExp(inSCE, altExpName), 
                                               X_name = useAssay,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               reducedDims = FALSE,
                                               metadata = FALSE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
    
  }
  
  if (method == "seurat") {
    scanpyObject$var['mean'] <- rowData(inSCE)$"Scanpy_mean"
    scanpyObject$var['std'] <- rowData(inSCE)$"Scanpy_std"
    sc$pp$highly_variable_genes(scanpyObject, 
                                flavor = method,
                                n_top_genes = as.integer(hvgNumber),
                                min_mean = minMean,
                                max_mean = maxMean,
                                min_disp = minDisp,
                                max_disp = maxDisp)
    
    #tmpSCE <- zellkonverter::AnnData2SCE(adata = scanpyObject)
    metadata(inSCE)$hvg <- scanpyObject$uns['hvg']
    rowData(inSCE)$dispersions <-
      unlist(scanpyObject$var['dispersions'])
    rowData(inSCE)$dispersions_norm <-
      unlist(scanpyObject$var['dispersions_norm'])
    rowData(inSCE)$means <-
      unlist(scanpyObject$var['means'])
    rowData(inSCE)$highly_variable <-
      unlist(scanpyObject$var['highly_variable'])
    
    metadata(inSCE)$sctk$runFeatureSelection$seurat <-
      list(
        useAssay = useAssay,
        rowData = c(
          "dispersions",
          "dispersions_norm",
          "means"
        )
      )
  } 
  #for this approach it is required that sce basic filering of cells and genes 
  #must be done
  else if (method == "cell_ranger") {
    sc$pp$highly_variable_genes(scanpyObject, 
                                flavor = method,
                                n_top_genes = as.integer(hvgNumber),
                                min_mean = minMean,
                                max_mean = maxMean,
                                min_disp = minDisp,
                                max_disp = maxDisp)
    
    
    #metadata(inSCE)$hvg <- metadata(tmpSCE)['hvg'][['hvg']] 
    metadata(fullSCE_cellranger)$hvg <- scanpyObject$uns['hvg']
    if (!altExp) {
      rowData(inSCE)$dispersions <-
        unlist(scanpyObject$var['dispersions'])
      rowData(inSCE)$dispersions_norm <-
        unlist(scanpyObject$var['dispersions_norm'])
      rowData(inSCE)$means <-
        unlist(scanpyObject$var['means'])
      rowData(inSCE)$highly_variable <-
        unlist(scanpyObject$var['highly_variable'])
      
    }
    else{
      altExpRows <- match(rownames(inSCE), rownames(scanpyObject))
      rowData(inSCE)$dispersions <-
        unlist(scanpyObject$var['dispersions'])[altExpRows]
      rowData(inSCE)$dispersions_norm <-
        unlist(scanpyObject$var['dispersions_norm'])[altExpRows]   
      rowData(inSCE)$means <-
        unlist(scanpyObject$var['means'])[altExpRows]
      rowData(inSCE)$highly_variable <-
        unlist(scanpyObject$var['highly_variable'])
    }
    
    metadata(inSCE)$sctk$runFeatureSelection$cell_ranger <-
      list(
        useAssay = useAssay,
        rowData = c(
          "dispersions",
          "dispersions_norm",
          "means",
          "highly_variable"
        )
      )
    
    mergedRowData <- merge(rowData(fullSCE_cellranger), rowData(inSCE)[, metadata(inSCE)$sctk$runFeatureSelection$cell_ranger$rowData],
                           by = 'row.names', all = TRUE)
    row.names(mergedRowData) <- mergedRowData$Row.names
    mergedRowData$Row.names <- NULL
    rowData(fullSCE_cellranger) <- mergedRowData
    inSCE <- fullSCE_cellranger
  }
  
  else if (method == "seurat_v3") {
    sc$pp$highly_variable_genes(scanpyObject, 
                                flavor = method,
                                n_top_genes = as.integer(hvgNumber),
                                min_mean = minMean,
                                max_mean = maxMean,
                                min_disp = minDisp,
                                max_disp = maxDisp)
    
    metadata(inSCE)$hvg <- scanpyObject$uns['hvg']
    
    rowData(inSCE)$variances <-
      unlist(scanpyObject$var['variances'])
    rowData(inSCE)$variances_norm <-
      unlist(scanpyObject$var['variances_norm'])   
    rowData(inSCE)$means <-
      unlist(scanpyObject$var['means'])
    rowData(inSCE)$highly_variable <-
      unlist(scanpyObject$var['highly_variable'])
    metadata(inSCE)$sctk$runFeatureSelection$seurat_v3 <-
      list(
        useAssay = useAssay,
        rowData = c(
          "variances",
          "variances_norm",
          "means"
        )
      )
  }
  return(inSCE)
}


#' plotScanpyHVG
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param log Plot on logarithmic axes. Default \code{FALSE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' plotScanpyHVG(sce)
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyHVG <- function(inSCE,
                          log = FALSE){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = FALSE,
                                             metadata = TRUE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = TRUE,
                                             verbose = NULL)
  if(is.null(scanpyObject$uns['hvg'])){
    stop(
      " High variable genes not found. Please run the 'runScanpyFindHVG' first."
    )
  }
  return(sc$pl$highly_variable_genes(scanpyObject,
                                     log = log))
}

################################################
######## PCA Function #########################
################################################

#' runScanpyPCA
#' Computes PCA on the input sce object and stores the calculated principal
#' components within the sce object
#' @param inSCE (sce) object on which to compute PCA
#' @param useAssay Assay containing scaled counts to use in PCA. Default
#' \code{"scanpyScaledData"}.
#' @param reducedDimName Name of new reducedDims object containing Scanpy PCA.
#' Default \code{scanpyPCA}.
#' @param nPCs numeric value of how many components to compute. Default
#' \code{50}.
#' @param method selected method to use for computation of pca. 
#' One of \code{'arpack'}, \code{'randomized'}, \code{'auto'} or \code{'lobpcg'}.
#' Default \code{"arpack"}.
#' @param use_highly_variable boolean value of whether to use highly variable 
#' genes only. By default uses them if they have been determined beforehand. 
#' @param seed Specify numeric value to set as a seed. Default \code{12345}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' computed principal components
#' @export
#' @importFrom SingleCellExperiment reducedDim<- rowSubset
#' @importFrom S4Vectors metadata<-
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyPCA <- function(inSCE,
                         useAssay = "scanpyScaledData",
                         reducedDimName = "scanpyPCA",
                         nPCs = 50,
                         method = c("arpack", "randomized", "auto", "lobpcg"),
                         use_highly_variable = TRUE,
                         seed = 12345){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  method <- match.arg(method)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = FALSE,
                                             metadata = FALSE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = FALSE,
                                             verbose = NULL)
  sc$tl$pca(scanpyObject, 
            svd_solver= method, 
            n_comps = as.integer(nPCs), 
            use_highly_variable = use_highly_variable)
  
  metadata(inSCE)$scanpy$PCA <- scanpyObject
  
  temp <- scanpyObject$obsm[['X_pca']]
  #colnames(temp) <- paste0(rep("PC", ncol(temp)), seq(ncol(10)))
  rownames(temp) <- colnames(inSCE)
  rotation <- scanpyObject$varm[['PCs']]
  rownames(rotation) <- rownames(inSCE)
  
  reducedDim(inSCE, reducedDimName) <- temp
  attr(reducedDim(inSCE, reducedDimName), "percentVar") <- 
    scanpyObject$uns[['pca']][['variance_ratio']]
  attr(reducedDim(inSCE, reducedDimName), "rotation") <- rotation
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  
  return(inSCE)
}


#' plotScanpyPCA
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param reducedDimName Name of new reducedDims object containing Scanpy PCA. 
#' @param color Keys for annotations of observations/cells or variables/genes.
#' @param title Provide title for panels either as string or list of strings
#' @param legend Location of legend, either 'on data', 'right margin' or a 
#' valid keyword for the loc parameter of Legend.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' plotScanpyPCA(sce)
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyPCA <- function(inSCE,
                          reducedDimName = "scanpyPCA", 
                          color = NULL,
                          title = '',
                          legend = 'right margin'){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(!reducedDimName %in% reducedDimNames(inSCE)){
    stop(
      "PCA results not found. Please run the 'runScanpyPCA' function first."
    )
  }
  reducedDim (inSCE, "X_pca") <- reducedDim(inSCE, reducedDimName)
  useAssay = metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]]$useAssay
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = TRUE,
                                             metadata = FALSE,
                                             colPairs = TRUE,
                                             rowPairs = TRUE,
                                             verbose = NULL)
  return(sc$pl$pca(scanpyObject,
                   color = color,
                   title = title, 
                   legend_loc = legend))
}


#' plotScanpyPCAGeneRanking
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param PC_comp For example, '1,2,3' means [1, 2, 3], first, second, 
#' third principal component.
#' @param includeLowest Whether to show the variables with both highest and 
#' lowest loadings. Default \code{TRUE}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' plotScanpyPCAGeneRanking(sce)
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyPCAGeneRanking <- function(inSCE, 
                                     PC_comp = "1,2,3",
                                     includeLowest = TRUE){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- metadata(inSCE)$scanpy$PCA
  return(sc$pl$pca_loadings(scanpyObject,
                            components = PC_comp,
                            include_lowest = includeLowest))
  
}

#' plotScanpyPCAVariance
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param nPCs Number of PCs to show. Default \code{50}.
#' @param log Plot on logarithmic scale. Default \code{FALSE}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' plotScanpyPCAVariance(sce)
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyPCAVariance <- function(inSCE,
                                  nPCs = 50,
                                  log = FALSE){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- metadata(inSCE)$scanpy$PCA
  return(sc$pl$pca_variance_ratio(scanpyObject,
                                  n_pcs = as.integer(nPCs),
                                  log = log))
  
}

########################################################
###### Clustering function ############################
#######################################################

#' runScanpyFindClusters
#' Computes the clusters from the input sce object and stores them back in sce
#' object
#' @param inSCE (sce) object from which clusters should be computed and stored
#' in
#' @param useAssay Assay containing scaled counts to use for clustering.
#' @param useReducedDim Reduction method to use for computing clusters. 
#' Default \code{"scanpyPCA"}.
#' @param nNeighbors The size of local neighborhood (in terms of number of 
#' neighboring data points) used for manifold approximation. Larger values 
#' result in more global views of the manifold, while smaller values result in 
#' more local data being preserved. Default \code{10}.
#' @param dims numeric value of how many components to use for computing
#' clusters. Default \code{40}.
#' @param method selected method to compute clusters. One of "louvain",
#' and "leiden". Default \code{louvain}.
#' @param colDataName Specify the name to give to this clustering result. 
#'  Default is \code{NULL} that will generate a meaningful name automatically.
#' @param resolution A parameter value controlling the coarseness of the 
#' clustering. Higher values lead to more clusters Default \code{1}.
#' @param niterations How many iterations of the Leiden clustering method to 
#' perform. Positive values above 2 define the total number of iterations to 
#' perform, -1 has the method run until it reaches its optimal clustering.
#' Default \code{-1}.
#' @param flavor Choose between to packages for computing the clustering.
#' Default \code{vtraag} 
#' @param use_weights Boolean. Use weights from knn graph. Default \code{FALSE}
#' @param cor_method correlation method to use. Options are ‘pearson’, 
#' ‘kendall’, and ‘spearman’. Default \code{pearson}.
#' @param inplace If True, adds dendrogram information to annData object, 
#' else this function returns the information. Default \code{TRUE}
#' @param externalReduction Pass DimReduce object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @param seed Specify numeric value to set as a seed. Default \code{12345}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' }
#' @return Updated sce object which now contains the computed clusters
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyFindClusters <- function(inSCE,
                                  useAssay = "scanpyScaledData",
                                  useReducedDim = "scanpyPCA",
                                  nNeighbors = 10,
                                  dims = 40,
                                  method = c("leiden", "louvain"),
                                  colDataName = NULL,
                                  resolution = 1,
                                  niterations = -1,
                                  flavor = 'vtraag',
                                  use_weights = FALSE,
                                  cor_method = 'pearson',
                                  inplace = TRUE,
                                  externalReduction = NULL,
                                  seed = 12345) {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!reticulate::py_module_available(module = "louvain")) {
    warning("Cannot find python module 'louvain', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, louvain can be installed on the local machine",
            "with pip (e.g. pip install louvain) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!reticulate::py_module_available(module = "leidenalg")) {
    warning("Cannot find python module 'leidenalg', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, leidenalg can be installed on the local machine",
            "with pip (e.g. pip install leidenalg) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }
  
  method <- match.arg(method)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,  
                                             X_name = useAssay,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = useReducedDim,
                                             metadata = FALSE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = FALSE,
                                             verbose = NULL)
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReducedDim <- "pca"
  } 
  
  if(is.null(colDataName)){
    colDataName = paste0("Scanpy", "_", method, "_", resolution) 
  }
  
  sc$pp$neighbors(scanpyObject, 
                  n_neighbors = as.integer(nNeighbors), 
                  n_pcs = as.integer(dims),
                  use_rep = useReducedDim)
  
  if (method == "louvain") {
    sc$tl$louvain(adata = scanpyObject,
                  key_added = colDataName,
                  flavor = flavor,
                  use_weights = use_weights)
  } else if (method == "leiden") {
    sc$tl$leiden(adata = scanpyObject,
                 key_added = colDataName,
                 n_iterations = as.integer(niterations))
  } 
  
  colData(inSCE)[[colDataName]] <-
    as.factor(unlist(scanpyObject$obs[colDataName]))
  S4Vectors::metadata(inSCE)$scanpy[method] <- colDataName
  
  return(inSCE)
}

###############################################
###### Embedding functions #####################
###############################################

#' runScanpyUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back
#' into the sce object
#' @param inSCE (sce) object on which to compute the UMAP
#' @param useAssay Specify name of assay to use. Default is \code{NULL}, so
#' \code{useReducedDim} param will be used instead.
#' @param useReducedDim Reduction to use for computing UMAP. 
#' Default is \code{"scanpyPCA"}.
#' @param reducedDimName Name of new reducedDims object containing Scanpy UMAP
#' Default \code{scanpyUMAP}.
#' @param dims Numerical value of how many reduction components to use for UMAP
#' computation. Default \code{40}.
#' @param minDist Sets the \code{"min_dist"} parameter to the underlying UMAP
#' call. Default \code{0.5}.
#' @param nNeighbors Sets the \code{"n_neighbors"} parameter to the underlying
#' UMAP call. Default \code{10}.
#' @param spread Sets the \code{"spread"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param alpha Sets the \code{"alpha"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param gamma Sets the \code{"gamma"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param externalReduction Pass DimReduce object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @param seed Specify numeric value to set as a seed. Default \code{12345}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' }
#' @return Updated sce object with UMAP computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyUMAP <- function(inSCE,
                          useAssay = NULL,
                          useReducedDim = "scanpyPCA",
                          reducedDimName = "scanpyUMAP",
                          dims = 40,  
                          minDist = 0.5,
                          nNeighbors = 10, 
                          spread = 1,
                          alpha=1.0, 
                          gamma=1.0, 
                          externalReduction = NULL,
                          seed = 12345) {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  if(!is.null(useAssay)){
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                               X_name = useAssay,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               reducedDims = useReducedDim,
                                               metadata = FALSE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
  } else{
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                               reducedDims = useReducedDim,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               metadata = FALSE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
  }
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReducedDim <- "pca"
  } 
  
  sc$pp$neighbors(scanpyObject, 
                  n_neighbors = as.integer(nNeighbors), 
                  n_pcs = as.integer(dims),
                  use_rep = useReducedDim)
  
  
  sc$tl$umap(scanpyObject, 
             n_components = 2L,
             min_dist = minDist,
             alpha = alpha, 
             gamma = gamma,
             spread = spread)
  
  
  temp <- scanpyObject$obsm[['X_umap']]
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  
  return(inSCE)
}



#' runScanpyTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back
#' into the sce object
#' @param inSCE (sce) object on which to compute the tSNE
#' @param useAssay Specify name of assay to use. Default is \code{NULL}, so
#' \code{useReducedDim} param will be used instead.
#' @param useReducedDim selected reduction method to use for computing tSNE.
#' Default \code{"scanpyPCA"}.
#' @param reducedDimName Name of new reducedDims object containing Scanpy tSNE
#' Default \code{scanpyTSNE}.
#' @param dims Number of reduction components to use for tSNE computation.
#' Default \code{40}.
#' @param perplexity Adjust the perplexity tuneable parameter for the underlying
#' tSNE call. Default \code{30}.
#' @param externalReduction Pass DimReduc object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @param seed Specify numeric value to set as a seed. Default \code{12345}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyTSNE(sce, useReducedDim = "scanpyPCA")
#' }
#' @return Updated sce object with tSNE computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyTSNE <- function(inSCE,
                          useAssay = NULL,
                          useReducedDim = "scanpyPCA", 
                          reducedDimName = "scanpyTSNE",
                          dims = 40,
                          perplexity = 30,
                          externalReduction = NULL,
                          seed = 12345){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  if(!is.null(useAssay)){
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                               X_name = useAssay,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               reducedDims = useReducedDim,
                                               metadata = FALSE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
  } else{
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                               reducedDims = useReducedDim,
                                               assays = FALSE,
                                               colData = TRUE,
                                               rowData = TRUE,
                                               varm = TRUE,
                                               metadata = FALSE,
                                               colPairs = FALSE,
                                               rowPairs = FALSE,
                                               skip_assays = FALSE,
                                               verbose = NULL)
  }
  
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReducedDim <- "pca"
  }
  
  sc$tl$tsne(scanpyObject, 
             n_pcs = as.integer(dims),
             use_rep = useReducedDim, 
             perplexity = as.integer(perplexity))
  
  
  
  temp <- scanpyObject$obsm[['X_tsne']]
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  
  return(inSCE)
}

#' plotScanpyEmbedding
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param reducedDimName Name of reducedDims object containing embeddings.
#' Eg. scanpyUMAP.
#' @param useAssay Specify name of assay to use. Default is \code{NULL},
#' which will use scaled assay by default. 
#' @param color Keys for annotations of observations/cells or variables/genes.
#' @param title Provide title for panels either as string or list of strings
#' @param legend Location of legend, either 'on data', 'right margin' or a 
#' valid keyword for the loc parameter of Legend.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' plotScanpyEmbedding(sce, reducedDimName = "scanpyUMAP", color = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyEmbedding <- function(inSCE,
                                reducedDimName,
                                useAssay = NULL,
                                color = NULL,
                                legend = 'right margin',
                                title = ''){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(!reducedDimName %in% reducedDimNames(inSCE)){
    stop(
      "Embedding result not found. Please run the 'runScanpyUMAP' or 
      'runScanpyTSNE' first."
    )
  }
  useAssay <- metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]]$useAssay
  if(is.null(useAssay)){
    useAssay <- "scanpyScaledData"
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay,
                                             assays = FALSE,
                                             colData = TRUE,
                                             rowData = TRUE,
                                             varm = TRUE,
                                             reducedDims = reducedDimName,
                                             metadata = FALSE,
                                             colPairs = FALSE,
                                             rowPairs = FALSE,
                                             skip_assays = FALSE,
                                             verbose = NULL)
  
  return(sc$pl$embedding(scanpyObject,
                         basis = reducedDimName,
                         color = color,
                         legend_loc = legend,
                         title = title))
  
}

###################################################
########## Marker Genes function ###################
###################################################

#' runScanpyFindMarkers
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param nGenes The number of genes that appear in the returned tables. 
#' Defaults to all genes.
#' @param useAssay Specify the name of the assay to use for computation
#'  of marker genes. It is recommended to use log normalized assay. 
#' @param colDataName colData to use as the key of the observations grouping to 
#' consider.
#' @param group1 Name of group1. Subset of groups, to which comparison shall be 
#' restricted, or 'all' (default), for all groups.
#' @param group2 Name of group2. If 'rest', compare each group to the union of 
#' the rest of the group. If a group identifier, compare with respect to this 
#' group. Default is 'rest'
#' @param test Test to use for DE. Default \code{"t-test"}.
#' @param corr_method p-value correction method. Used only for 't-test', 
#' 't-test_overestim_var', and 'wilcoxon'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' }
#' @return A \code{SingleCellExperiment} object that contains marker genes
#' populated in a data.frame stored inside metadata slot.
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
runScanpyFindMarkers <- function(inSCE,
                                 nGenes = NULL,
                                 useAssay = "scanpyNormData",
                                 colDataName,
                                 group1 = "all",
                                 group2 = "rest",
                                 test = c("wilcoxon", "t-test", "t-test_overestim_var", "logreg"),
                                 corr_method = c("benjamini-hochberg", "bonferroni")) {
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  test <- match.arg(test)
  corr_method <- match.arg(corr_method)
  
  #store results in a temporary sce object
  tmpSCE <- SingleCellExperiment(
    assays = list(counts = counts(inSCE), temp = as.matrix(assay(inSCE, useAssay))))
  # store back colData from sce into the tmpSCE colData slot
  colData(tmpSCE)[colDataName] <- as.character(colData(inSCE)[[colDataName]])
  rowData(tmpSCE)$id <- rownames(inSCE)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = tmpSCE, X_name = "temp")
  if(!is.null(nGenes)){
    nGenes = as.integer(nGenes)
  }
  sc$tl$rank_genes_groups(scanpyObject, 
                          groupby = colDataName, 
                          groups = group1,
                          reference = group2,
                          method = test, 
                          n_genes = nGenes,
                          corr_method = corr_method)
  
  py <- reticulate::py
  py$scanpyObject <- scanpyObject
  reticulate::py_run_string("import pandas as pd")
  reticulate::py_run_string(
    "names = pd.DataFrame(scanpyObject.uns['rank_genes_groups']['names'])", 
    convert = TRUE)
  reticulate::py_run_string(
    "logFoldChanges = pd.DataFrame(scanpyObject.uns['rank_genes_groups']['logfoldchanges'])", 
    convert = TRUE)
  reticulate::py_run_string(
    "pvals_adj = pd.DataFrame(scanpyObject.uns['rank_genes_groups']['pvals_adj'])", 
    convert = TRUE)
  reticulate::py_run_string(
    "scores = pd.DataFrame(scanpyObject.uns['rank_genes_groups']['scores'])",
    convert = TRUE)
  
  markerGenesNames <- utils::stack(py$names)
  colnames(markerGenesNames) <- c("Gene","findMarker_cluster")
  Log2_FC <- unlist(py$logFoldChanges)
  Pvalue <- unlist(py$pvals_adj)
  zscore <- unlist(py$scores)
  
  
  markerGenesTable <- data.frame()
  markerGenesTable <- cbind(markerGenesNames, Log2_FC, Pvalue, zscore)
  
  S4Vectors::metadata(inSCE)$"findMarkerScanpyObject" <- scanpyObject
  S4Vectors::metadata(inSCE)$scanpyMarkersTable <- markerGenesTable
  return(inSCE)
}

#' plotScanpyMarkerGenes
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param groups The groups for which to show the gene ranking. Default \code{NULL}
#' means that all groups will be considered. 
#' @param nGenes Number of genes to show. Default \code{10}
#' @param nCols Number of panels shown per row. Default \code{4}
#' @param sharey Controls if the y-axis of each panels should be shared. 
#' Default \code{FALSE} allows each panel to have its own y-axis range.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenes(sce, groups = '0')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMarkerGenes <- function(inSCE,
                                  groups = NULL,
                                  nGenes = 10,
                                  nCols = 4,
                                  sharey = FALSE){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(is.null(metadata(inSCE)[["findMarkerScanpyObject"]])){
    stop("Marker genes not found. Please run the 'runScanpyFindMarkers' function first.")
  }
  scanpyObject <- metadata(inSCE)[["findMarkerScanpyObject"]]
  return(sc$pl$rank_genes_groups(scanpyObject, 
                                 groups = groups, 
                                 n_genes = as.integer(nGenes),
                                 ncols = as.integer(nCols), 
                                 sharey = sharey))
  
}

#' plotScanpyMarkerGenesViolin
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param groups The groups for which to show the gene ranking. Default \code{NULL}
#' means that all groups will be considered. 
#' @param features List of genes to plot. Is only useful if interested in a 
#' custom gene list
#' @param nGenes Number of genes to show. Default \code{10}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesViolin(sce, groups = '0')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMarkerGenesViolin <- function(inSCE,
                                        groups = NULL,
                                        features = NULL,
                                        nGenes = 10){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(is.null(metadata(inSCE)["findMarkerScanpyObject"])){
    stop("Marker genes not found. Please run the 'runScanpyFindMarkers' function first.")
  }
  scanpyObject <- metadata(inSCE)[["findMarkerScanpyObject"]]
  return(sc$pl$rank_genes_groups_violin(scanpyObject,
                                        groups = groups,
                                        gene_names = features,
                                        n_genes = as.integer(nGenes)))
  
}

#' plotScanpyMarkerGenesHeatmap
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param groups The groups for which to show the gene ranking. Default \code{NULL}
#' means that all groups will be considered.
#' @param groupBy The key of the observation grouping to consider. By default, 
#' the groupby is chosen from the rank genes groups parameter.
#' @param nGenes Number of genes to show. Default \code{10}
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes). The var_names could be a dictionary or a list. 
#' @param log2fcThreshold Only output DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesHeatmap(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMarkerGenesHeatmap <- function(inSCE,
                                         groups = NULL,
                                         groupBy,
                                         nGenes = 10,
                                         features = NULL,
                                         log2fcThreshold = NULL){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(is.null(metadata(inSCE)["findMarkerScanpyObject"])){
    stop("Marker genes not found. Please run the 'runScanpyFindMarkers' function first.")
  }
  scanpyObject <- metadata(inSCE)[["findMarkerScanpyObject"]]
  
  if(!is.null(features)){
    return(sc$pl$rank_genes_groups_heatmap(scanpyObject,
                                           groups = groups,
                                           groupby = groupBy,
                                           var_names = features,
                                           min_logfoldchange = log2fcThreshold,
                                           show_gene_labels = TRUE,
                                           dendrogram = FALSE))
  }
  else
    return(sc$pl$rank_genes_groups_heatmap(scanpyObject,
                                           groups = groups,
                                           groupby = groupBy,
                                           n_genes = as.integer(nGenes),
                                           var_names = NULL,
                                           min_logfoldchange = log2fcThreshold,
                                           show_gene_labels = TRUE,
                                           dendrogram = FALSE))
  
}

#' plotScanpyMarkerGenesDotPlot
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param groups The groups for which to show the gene ranking. Default \code{NULL}
#' means that all groups will be considered.
#' @param nGenes Number of genes to show. Default \code{10}
#' @param groupBy The key of the observation grouping to consider. By default, 
#' the groupby is chosen from the rank genes groups parameter.
#' @param log2fcThreshold Only output DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}.
#' @param parameters The options for marker genes results to plot are: 
#' ‘scores’, ‘logfoldchanges’, ‘pvals’, ‘pvals_adj’, ‘log10_pvals’, ‘log10_pvals_adj’.
#' If NULL provided then it uses mean gene value to plot.  
#' @param standardScale Whether or not to standardize the given dimension 
#' between 0 and 1, meaning for each variable or group, subtract the minimum and 
#' divide each by its maximum. Default \code{NULL} means that it doesn't perform
#' any scaling. 
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes) to check their fold changes or p-values, instead of 
#'  the top/bottom genes. The gene names could be a dictionary or a list. 
#'  Default \code{NULL}
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin. 
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param colorBarTitle Title for the color bar. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesDotPlot(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMarkerGenesDotPlot <- function(inSCE,
                                         groups = NULL,
                                         nGenes = 10, 
                                         groupBy,
                                         log2fcThreshold = NULL,
                                         parameters = "logfoldchanges",
                                         standardScale = NULL,
                                         features = NULL,
                                         title = '',
                                         vmin = NULL,
                                         vmax = NULL,
                                         colorBarTitle = "log fold change"){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(is.null(metadata(inSCE)["findMarkerScanpyObject"])){
    stop("Marker genes not found. Please run the 'runScanpyFindMarkers' function first.")
  }
  scanpyObject <- metadata(inSCE)[["findMarkerScanpyObject"]]
  
  if(!is.null(features)){
    return(sc$pl$rank_genes_groups_dotplot(scanpyObject,
                                           groups = groups,
                                           groupby = groupBy,
                                           min_logfoldchange = log2fcThreshold,
                                           values_to_plot = parameters,
                                           standard_scale = standardScale,
                                           var_names = features,
                                           title = title,
                                           vmin = vmin,
                                           vmax = vmax,
                                           cmap = 'bwr',
                                           dendrogram = FALSE,
                                           colorbar_title = colorBarTitle))
  }
  else
    return(sc$pl$rank_genes_groups_dotplot(scanpyObject,
                                           groups = groups,
                                           n_genes = as.integer(nGenes),
                                           groupby = groupBy,
                                           min_logfoldchange = log2fcThreshold,
                                           values_to_plot = parameters,
                                           standard_scale = standardScale,
                                           var_names = NULL,
                                           title = title,
                                           vmin = vmin,
                                           vmax = vmax,
                                           cmap = 'bwr',
                                           dendrogram = FALSE,
                                           colorbar_title = colorBarTitle))
  
  
  
}

#' plotScanpyMarkerGenesMatrixPlot
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param groups The groups for which to show the gene ranking. Default \code{NULL}
#' means that all groups will be considered.
#' @param nGenes Number of genes to show. Default \code{10}
#' @param groupBy The key of the observation grouping to consider. By default, 
#' the groupby is chosen from the rank genes groups parameter.
#' @param log2fcThreshold Only output DEGs with the absolute values of log2FC
#' larger than this value. Default \code{NULL}.
#' @param parameters The options for marker genes results to plot are: 
#' ‘scores’, ‘logfoldchanges’, ‘pvals’, ‘pvals_adj’, ‘log10_pvals’, ‘log10_pvals_adj’.
#' If NULL provided then it uses mean gene value to plot.  
#' @param standardScale Whether or not to standardize the given dimension 
#' between 0 and 1, meaning for each variable or group, subtract the minimum and 
#' divide each by its maximum. Default \code{NULL} means that it doesn't perform
#' any scaling. 
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes) to check their fold changes or p-values, instead of 
#'  the top/bottom genes. The var_names could be a dictionary or a list. 
#'  Default \code{NULL}
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin. 
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param colorBarTitle Title for the color bar. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesMatrixPlot(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMarkerGenesMatrixPlot <- function(inSCE,
                                            groups = NULL,
                                            nGenes = 10, 
                                            groupBy,
                                            log2fcThreshold = NULL,
                                            parameters = "logfoldchanges",
                                            standardScale = 'var',
                                            features = NULL,
                                            title = '',
                                            vmin = NULL,
                                            vmax = NULL,
                                            colorBarTitle = "log fold change"){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if(is.null(metadata(inSCE)["findMarkerScanpyObject"])){
    stop("Marker genes not found. Please run the 'runScanpyFindMarkers' function first.")
  }
  scanpyObject <- metadata(inSCE)[["findMarkerScanpyObject"]]
  
  if(!is.null(features)){
    return(sc$pl$rank_genes_groups_matrixplot(scanpyObject,
                                              groups = groups,
                                              groupby = groupBy,
                                              min_logfoldchange = log2fcThreshold,
                                              values_to_plot = parameters,
                                              standard_scale = standardScale,
                                              var_names = features,
                                              title = title,
                                              vmin = vmin,
                                              vmax = vmax,
                                              cmap = 'bwr',
                                              dendrogram = FALSE,
                                              colorbar_title = colorBarTitle))
  }
  else
    return(sc$pl$rank_genes_groups_matrixplot(scanpyObject,
                                              groups = groups,
                                              n_genes = as.integer(nGenes),
                                              groupby = groupBy,
                                              min_logfoldchange = log2fcThreshold,
                                              values_to_plot = parameters,
                                              standard_scale = standardScale,
                                              var_names = NULL,
                                              title = title,
                                              vmin = vmin,
                                              vmax = vmax,
                                              cmap = 'bwr',
                                              dendrogram = FALSE,
                                              colorbar_title = colorBarTitle))
  
  
}



###############################################
####### General plotting functions ############
################################################

#' plotScanpyHeatmap
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Assay to use for plotting. By default it will use counts
#' assay.
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes). The var_names could be a dictionary or a list. 
#' @param groupBy The key of the observation grouping to consider.
#' @param standardScale Whether or not to standardize the given dimension 
#' between 0 and 1, meaning for each variable or group, subtract the minimum and 
#' divide each by its maximum. Default \code{NULL} means that it doesn't perform
#' any scaling. 
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin.
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyHeatmap(sce, features = markers, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyHeatmap <- function(inSCE,
                              useAssay = NULL,
                              features,
                              groupBy,
                              standardScale = 'var',
                              vmin = NULL,
                              vmax = NULL){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay)
  
  return(sc$pl$heatmap(scanpyObject,
                       var_names = features,
                       groupby = groupBy,
                       standard_scale = standardScale,
                       vmin = vmin,
                       vmax = vmax,
                       show_gene_labels = TRUE,
                       dendrogram = FALSE))
  
}


#' plotScanpyDotPlot
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Assay to use for plotting. By default it will use counts assay.
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes). The var_names could be a dictionary or a list. 
#' @param groupBy The key of the observation grouping to consider.
#' @param standardScale Whether or not to standardize the given dimension 
#' between 0 and 1, meaning for each variable or group, subtract the minimum and 
#' divide each by its maximum. Default \code{NULL} means that it doesn't perform
#' any scaling. 
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin.
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin. 
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param colorBarTitle Title for the color bar. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyDotPlot(sce, features = markers, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyDotPlot <- function(inSCE,
                              useAssay = NULL,
                              features,
                              groupBy,
                              standardScale = NULL,
                              title = '',
                              vmin = NULL,
                              vmax = NULL,
                              colorBarTitle = "Mean expression in group"){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay)
  
  return(sc$pl$dotplot(scanpyObject,
                       var_names = features,
                       groupby = groupBy,
                       standard_scale = standardScale,
                       title = title,
                       vmin = vmin,
                       vmax = vmax,
                       cmap = 'bwr',
                       dendrogram = FALSE,
                       colorbar_title = colorBarTitle))
  
  
}

#' plotScanpyMatrixPlot
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Assay to use for plotting. By default it will use counts assay.
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes). The var_names could be a dictionary or a list. 
#' @param groupBy The key of the observation grouping to consider.
#' @param standardScale Whether or not to standardize the given dimension 
#' between 0 and 1, meaning for each variable or group, subtract the minimum and 
#' divide each by its maximum. Default \code{NULL} means that it doesn't perform
#' any scaling. 
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin.
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param title Provide title for the figure.
#' @param vmin The value representing the lower limit of the color scale. 
#' Values smaller than vmin are plotted with the same color as vmin. 
#' Default \code{NULL}
#' @param vmax The value representing the upper limit of the color scale. 
#' Values larger than vmax are plotted with the same color as vmax. 
#' Default \code{NULL}
#' @param colorBarTitle Title for the color bar. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyMatrixPlot(sce, features = markers, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyMatrixPlot <- function(inSCE,
                              useAssay = NULL,
                              features,
                              groupBy,
                              standardScale = NULL,
                              title = '',
                              vmin = NULL,
                              vmax = NULL,
                              colorBarTitle = "Mean expression in group"){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay)
  
  return(sc$pl$matrixplot(scanpyObject,
                       var_names = features,
                       groupby = groupBy,
                       standard_scale = standardScale,
                       title = title,
                       vmin = vmin,
                       vmax = vmax,
                       cmap = 'bwr',
                       dendrogram = FALSE,
                       colorbar_title = colorBarTitle))
  
  
}


#' plotScanpyViolin
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Assay to use for plotting. By default it will use counts
#' assay.
#' @param features Genes to plot. Sometimes is useful to pass a specific list of
#'  var names (e.g. genes). The var_names could be a dictionary or a list. 
#' @param groupBy The key of the observation grouping to consider.
#' @param xlabel Label of the x axis. Defaults to groupBy.
#' @param ylabel Label of the y axis. If NULL and groupBy is NULL, 
#' defaults to 'value'. If NULL and groupBy is not NULL, defaults to features.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyNormData", method = "seurat")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyScaledData")
#' sce <- runScanpyFindClusters(sce, useReducedDim = "scanpyPCA")
#' sce <- runScanpyUMAP(sce, useReducedDim = "scanpyPCA")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyViolin(sce, features = markers, groupBy = "Scanpy_louvain_1")
#' }
#' @return plot object
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
plotScanpyViolin <- function(inSCE,
                             useAssay = NULL,
                             features,
                             groupBy, 
                             xlabel = '', 
                             ylabel = NULL){
  
  if (!reticulate::py_module_available(module = "scanpy")) {
    warning("Cannot find python module 'scanpy', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scanpy can be installed on the local machine",
            "with pip (e.g. pip install scanpy) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,
                                             X_name = useAssay)
  
  return(sc$pl$violin(scanpyObject,
                      keys = features,
                      groupby = groupBy,
                      xlabel = xlabel,
                      ylabel = ylabel))
  
}
