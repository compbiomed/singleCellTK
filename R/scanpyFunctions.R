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
#' before normalization. Default \code{1}
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
runScanpyNormalizeData <- function(inSCE,
                                   useAssay,
                                   targetSum = 1,
                                   maxFraction = 0.05,
                                   normAssayName = "scanpyNormData") {
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  normValue <- sc$pp$normalize_total(scanpyObject, 
                                     target_sum = targetSum, 
                                     max_fraction =  maxFraction,
                                     inplace = FALSE)
  
  scanpyObject$obs$n_counts <- normValue$norm_factor
  scanpyObject$X <- normValue$X
  
  # log transform the data.
  sc$pp$log1p(scanpyObject, copy = TRUE)
  
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, normAssayName)
  metadata(inSCE)$scanpy$obj <- scanpyObject
  metadata(inSCE)$scanpy$normValues <- normValue$X
  inSCE@metadata$scanpy$normAssay <- normAssayName
  colData(inSCE)$n_counts <- 
    as.factor(unlist(metadata(inSCE)$scanpy$obj$obs['n_counts']))
  
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
runScanpyScaleData <- function(inSCE,
                               useAssay = "scanpyNormData",
                               scaledAssayName = "scanpyScaledData") {
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  sc$pp$scale(scanpyObject, max_value=10)
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, scaledAssayName, "X")
  metadata(inSCE)$scanpy$obj <- scanpyObject
  metadata(inSCE)$scanpy$scaledAssay <- scaledAssayName
  rowData(inSCE)$"Scanpy_mean" <- metadata(inSCE)$scanpy$obj$var$mean 
  rowData(inSCE)$"Scanpy_std" <- metadata(inSCE)$scanpy$obj$var$std 
  
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
#'  of variable genes. It is recommended to use a raw counts assay with the
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
#' flavor='seurat_v3'.
#' @param maxMean If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'.
#' @param minDisp If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'.
#' @param maxDisp If n_top_genes unequals None, this and all other cutoffs for 
#' the means and the normalized dispersions are ignored. Ignored if 
#' flavor='seurat_v3'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindHVG(sce, useAssay = "scanpyScaledData", method = "seurat")
#' g <- getTopHVG(sce, method = "seurat", hvgNumber = 500)
#' }
#' @return Updated \code{SingleCellExperiment} object with highly variable genes
#' computation stored
#' \code{\link{getTopHVG}}, \code{\link{plotTopHVG}}
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors metadata<-
runScanpyFindHVG <- function(inSCE,
                             useAssay = "scanpyScaledData",
                             method = c("seurat", "cell_ranger", "seurat_v3"),
                             altExpName = "featureSubset",
                             altExp = FALSE,
                             hvgNumber = 2000,
                             minMean = 0.0125,
                             maxMean = 3,
                             minDisp = 0.5,
                             maxDisp = Inf) {
  
  method <- match.arg(method)
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  if (!altExp) {
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, 
                                               X_name = useAssay)
  }
  else{
    scanpyObject <- zellkonverter::SCE2AnnData(sce = altExp(inSCE, altExpName), 
                                               X_name = useAssay)
    
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
    
    metadata(inSCE)$scanpy$hvg <- scanpyObject
    rowData(inSCE)$scanpy_variableFeatures_seurat_dispersion <-
      inSCE@metadata$scanpy$hvg["var"][["dispersions"]]
    rowData(inSCE)$scanpy_variableFeatures_seurat_dispersionScaled <-
      inSCE@metadata$scanpy$hvg["var"][["dispersions_norm"]]    
    rowData(inSCE)$scanpy_variableFeatures_seurat_mean <-
      inSCE@metadata$scanpy$hvg["var"][["means"]]  
    
    metadata(inSCE)$sctk$runFeatureSelection$seurat <-
      list(
        useAssay = useAssay,
        rowData = c(
          "scanpy_variableFeatures_seurat_dispersion",
          "scanpy_variableFeatures_seurat_dispersionScaled",
          "scanpy_variableFeatures_seurat_mean"
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
    
    metadata(inSCE)$scanpy$hvg <- scanpyObject
    if (!altExp) {
      rowData(inSCE)$scanpy_variableFeatures_cr_dispersion <-
        inSCE@metadata$scanpy$hvg["var"][["dispersions"]]
      rowData(inSCE)$scanpy_variableFeatures_cr_dispersionScaled <-
        inSCE@metadata$scanpy$hvg["var"][["dispersions_norm"]]    
      rowData(inSCE)$scanpy_variableFeatures_cr_mean <-
        inSCE@metadata$scanpy$hvg["var"][["means"]] 
      
    }
    else{
      scanpyToSCE <- zellkonverter::AnnData2SCE(adata = scanpyObject)
      altExpRows <- match(rownames(inSCE), rownames(scanpyToSCE))
      rowData(inSCE)$scanpy_variableFeatures_cr_dispersion <-
        inSCE@metadata$scanpy$hvg["var"][["dispersions"]][altExpRows]
      rowData(inSCE)$scanpy_variableFeatures_cr_dispersionScaled <-
        inSCE@metadata$scanpy$hvg["var"][["dispersions_norm"]] [altExpRows]   
      rowData(inSCE)$scanpy_variableFeatures_cr_mean <-
        inSCE@metadata$scanpy$hvg["var"][["means"]][altExpRows]
    }
    
    metadata(inSCE)$sctk$runFeatureSelection$cell_ranger <-
      list(
        useAssay = useAssay,
        rowData = c(
          "scanpy_variableFeatures_cr_dispersion",
          "scanpy_variableFeatures_cr_dispersionScaled",
          "scanpy_variableFeatures_cr_mean"
        )
      )
  }
  
  else if (method == "seurat_v3") {
    sc$pp$highly_variable_genes(scanpyObject, 
                                flavor = method,
                                n_top_genes = as.integer(hvgNumber),
                                min_mean = minMean,
                                max_mean = maxMean,
                                min_disp = minDisp,
                                max_disp = maxDisp)
    
    metadata(inSCE)$scanpy$hvg <- scanpyObject
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_variances <-
      inSCE@metadata$scanpy$hvg["var"][["variances"]]
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_variancesScaled <-
      inSCE@metadata$scanpy$hvg["var"][["variances_norm"]]    
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_mean <-
      inSCE@metadata$scanpy$hvg["var"][["means"]]  
    metadata(inSCE)$sctk$runFeatureSelection$seurat_v3 <-
      list(
        useAssay = useAssay,
        rowData = c(
          "scanpy_variableFeatures_seuratv3_variances",
          "scanpy_variableFeatures_seuratv3_variancesScaled",
          "scanpy_variableFeatures_seuratv3_mean"
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
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' plotScanpyHVG(sce)
#' }
#' @return plot object
#' @export
#findhvg stored in scanpy$hvg rather than scanpy$obj. plot function made
plotScanpyHVG <- function(inSCE,
                          log = FALSE){
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  if(is.null(scanpyObject$uns['scanpy']$hvg)){
    stop(
      " High variable genes not found. Kindly run 'runScanpyFindHVG' first "
    )
  }
  return(sc$pl$highly_variable_genes(scanpyObject$uns['scanpy']$hvg,
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
#' \code{20}.
#' @param algorithm selected method to use for computation of pca. 
#' One of \code{'arpack'}, \code{'randomized'}, \code{'auto'} or \code{'lobpcg'}.
#' Default \code{"arpack"}.
#' @param use_highly_variable boolean value of whether to use highly variable 
#' genes only. By default uses them if they have been determined beforehand. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' computed principal components
#' @export
#' @importFrom SingleCellExperiment reducedDim<- rowSubset
#' @importFrom S4Vectors metadata<-
#' 
runScanpyPCA <- function(inSCE,
                         useAssay = "scanpyNormData",
                         reducedDimName = "scanpyPCA",
                         nPCs = 50L,
                         algorithm = c("arpack", "randomized", "auto", "lobpcg"),
                         use_highly_variable = TRUE){
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  algorithm <- match.arg(algorithm)
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  if(use_highly_variable == TRUE){
    scanpyObject$var['highly_variable'] <- 
      scanpyObject$uns['scanpy']$hvg$var$highly_variable
  }
  
  sc$tl$pca(scanpyObject, 
            svd_solver= algorithm, 
            n_comps = nPCs, 
            use_highly_variable = use_highly_variable)
  
  metadata(inSCE)$scanpy$PCA <- scanpyObject
  
  temp <- scanpyObject$obsm['X_pca']
  #colnames(temp) <- paste0(rep("PC", ncol(temp)), seq(ncol(10)))
  rownames(temp) <- colnames(inSCE)
  rotation <- scanpyObject$varm['PCs']
  rownames(rotation) <- rownames(inSCE)
  
  reducedDim(inSCE, reducedDimName) <- temp
  attr(reducedDim(inSCE, reducedDimName), "percentVar") <- 
    scanpyObject$uns['pca'][['variance_ratio']]
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' plotScanpyPCA(sce)
#' }
#' @return plot object
#' @export
plotScanpyPCA <- function(inSCE,
                          reducedDimName = "scanpyPCA", 
                          color = NULL,
                          title = '',
                          legend = 'right margin'){
  
  if(!reducedDimName %in% reducedDimNames(inSCE)){
    stop(
      "PCA results not found. Perform runScanpyPCA first. "
    )
  }
  #As scanpyPCA function will find the X_pca variable where it originally 
  #stores the result, thereby we should have the results saved by that variable
  #name
  reducedDim (inSCE, "X_pca") <- reducedDim(inSCE, reducedDimName)
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' plotScanpyPCAGeneRanking(sce)
#' }
#' @return plot object
#' @export
plotScanpyPCAGeneRanking <- function(inSCE, 
                                     PC_comp = "1,2,3",
                                     includeLowest = TRUE){
  
  scanpyObject <- metadata(inSCE)$scanpy$PCA
  return(sc$pl$pca_loadings(scanpyObject,
                            components = PC_comp,
                            include_lowest = includeLowest))
  
}

#' plotScanpyPCAVariance
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param nPCs Number of PCs to show. Default \code{20}.
#' @param log Plot on logarithmic scale. Default \code{FALSE}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' plotScanpyPCAVariance(sce)
#' }
#' @return plot object
#' @export
plotScanpyPCAVariance <- function(inSCE,
                                  nPCs = 20,
                                  log = FALSE){
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
#' @param useReduction Reduction method to use for computing clusters. 
#' Default \code{"pca"}.
#' @param nNeighbors The size of local neighborhood (in terms of number of 
#' neighboring data points) used for manifold approximation. Larger values 
#' result in more global views of the manifold, while smaller values result in 
#' more local data being preserved. Default \code{15L}.
#' @param dims numeric value of how many components to use for computing
#' clusters. Default \code{10}.
#' @param algorithm selected algorithm to compute clusters. One of "louvain",
#' and "leiden". Default \code{louvain}.
#' @param colDataName Specify the name to give to this clustering result. 
#'  Default is \code{NULL} that will generate a meaningful name automatically.
#' @param resolution A parameter value controlling the coarseness of the 
#' clustering. Higher values lead to more clusters Default \code{1}.
#' @param niterations How many iterations of the Leiden clustering algorithm to 
#' perform. Positive values above 2 define the total number of iterations to 
#' perform, -1 has the algorithm run until it reaches its optimal clustering.
#' @param flavor Choose between to packages for computing the clustering. 
#' @param use_weights Boolean. Use weights from knn graph.
#' @param cor_method correlation method to use. Options are ‘pearson’, 
#' ‘kendall’, and ‘spearman’. Default 'pearson'.
#' @param inplace If True, adds dendrogram information to annData object, 
#' else this function returns the information.
#' @param externalReduction Pass DimReduce object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' }
#' @return Updated sce object which now contains the computed clusters
#' @export
runScanpyFindClusters <- function(inSCE,
                                  useAssay = "scanpyScaledData",
                                  useReduction = "scanpyPCA",
                                  nNeighbors = 15L,
                                  dims = 2L,
                                  algorithm = c("louvain", "leiden"),
                                  colDataName = NULL,
                                  resolution = 1,
                                  niterations = -1,
                                  flavor = 'vtraag',
                                  use_weights = FALSE,
                                  cor_method = 'pearson',
                                  inplace = TRUE,
                                  externalReduction = NULL) {
  algorithm <- match.arg(algorithm)
  useReduction <- useReduction
  if (missing(useAssay)) {
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message(
      "'useAssay' parameter missing. Using the first available assay ",
      "instead: '",
      useAssay,
      "'"
    )
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE,  X_name = useAssay)
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReduction <- "pca"
  } 
  
  if(is.null(colDataName)){
    colDataName = paste0("Scanpy", "_", algorithm, "_", resolution) 
  }
  
  sc$pp$neighbors(scanpyObject, 
                  n_neighbors = nNeighbors, 
                  n_pcs = dims,
                  use_rep = useReduction)
  
  if (algorithm == "louvain") {
    sc$tl$louvain(adata = scanpyObject,
                  key_added = colDataName,
                  flavor = flavor,
                  use_weights = use_weights)
  } else if (algorithm == "leiden") {
    sc$tl$leiden(adata = scanpyObject,
                 key_added = colDataName,
                 n_iterations = niterations)
  } 
  
  metadata(inSCE)$scanpy$obj <- scanpyObject
  colData(inSCE)[[colDataName]] <-
    as.factor(unlist(scanpyObject$obs[colDataName]))
  S4Vectors::metadata(inSCE)$scanpy[algorithm] <- colDataName
  
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
#' \code{useReduction} param will be used instead.
#' @param useReduction Reduction to use for computing UMAP. 
#' Default is \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Scanpy UMAP
#' Default \code{scanpyUMAP}.
#' @param dims Numerical value of how many reduction components to use for UMAP
#' computation. Default \code{10}.
#' @param minDist Sets the \code{"min_dist"} parameter to the underlying UMAP
#' call. Default \code{0.5}.
#' @param nNeighbors Sets the \code{"n_neighbors"} parameter to the underlying
#' UMAP call. Default \code{15L}.
#' @param spread Sets the \code{"spread"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param alpha Sets the \code{"alpha"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param gamma Sets the \code{"gamma"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' @param externalReduction Pass DimReduce object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' }
#' @return Updated sce object with UMAP computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
runScanpyUMAP <- function(inSCE,
                          useAssay = NULL,
                          useReduction = "scanpyPCA",
                          reducedDimName = "scanpyUMAP",
                          dims = 10L,  
                          minDist = 0.5,
                          nNeighbors = 15L, 
                          spread = 1,
                          alpha=1.0, 
                          gamma=1.0, 
                          externalReduction = NULL) {
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  if(!is.null(useAssay)){
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  } else{
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, reducedDims = useReduction)
  }
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReduction <- "pca"
  } 
  
  if(!is.null(useAssay)){
    sc$pp$neighbors(scanpyObject, 
                    n_neighbors = nNeighbors, 
                    n_pcs = 0)
  } else{
    sc$pp$neighbors(scanpyObject, 
                    n_neighbors = nNeighbors, 
                    n_pcs = dims,
                    use_rep = useReduction)
  }
  
  sc$tl$umap(scanpyObject, 
             n_components = dims,
             min_dist = minDist,
             alpha = alpha, 
             gamma = gamma,
             spread = spread)
  
  metadata(inSCE)$scanpy$obj <- scanpyObject
  
  temp <- scanpyObject$obsm['X_umap']
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
#' \code{useReduction} param will be used instead.
#' @param useReduction selected reduction algorithm to use for computing tSNE.
#' Default \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Scanpy tSNE
#' Default \code{scanpyTSNE}.
#' @param dims Number of reduction components to use for tSNE computation.
#' Default \code{10}.
#' @param perplexity Adjust the perplexity tuneable parameter for the underlying
#' tSNE call. Default \code{15}.
#' @param externalReduction Pass DimReduc object if PCA computed through
#' other libraries. Default \code{NULL}.
#' @return Updated sce object with tSNE computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
runScanpyTSNE <- function(inSCE,
                          useAssay = NULL,
                          useReduction = "scanpyPCA", 
                          reducedDimName = "scanpyTSNE",
                          dims = 10L,
                          perplexity = 15L,
                          externalReduction = NULL){
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  if(!is.null(useAssay)){
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  } else{
    scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, reducedDims = useReduction)
  }
  
  
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReduction <- "pca"
  }
  
  if(!is.null(useAssay)){
    sc$tl$tsne(scanpyObject, 
               n_pcs = dims,
               use_rep = 'X',
               perplexity = perplexity)
  }
  else{
    sc$tl$tsne(scanpyObject, 
               n_pcs = dims,
               use_rep = useReduction, 
               perplexity = perplexity)
  }
  
  metadata(inSCE)$scanpy$obj <- scanpyObject
  
  temp <- scanpyObject$obsm['X_tsne']
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
#' @param color Keys for annotations of observations/cells or variables/genes.
#' @param title Provide title for panels either as string or list of strings
#' @param legend Location of legend, either 'on data', 'right margin' or a 
#' valid keyword for the loc parameter of Legend.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' plotScanpyEmbedding(sce, reducedDimName = "scanpyUMAP", color = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
plotScanpyEmbedding <- function(inSCE,
                                reducedDimName,
                                color = NULL,
                                legend = 'right margin',
                                title = ''){
  
  if(!reducedDimName %in% reducedDimNames(inSCE)){
    stop(
      "Embedding result not found. Kindly implement embedding algorithms first"
    )
  }
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  
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
#'  of marker genes. It is recommended to use scaled assay. 
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
#' sce <- runScanpyScaleData(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts", algorithm = "louvain")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' }
#' @return A \code{SingleCellExperiment} object that contains marker genes
#' populated in a data.frame stored inside metadata slot.
#' @export
runScanpyFindMarkers <- function(inSCE,
                                 nGenes = NULL,
                                 useAssay = "scanpyScaledData",
                                 colDataName,
                                 group1 = "all",
                                 group2 = "rest",
                                 test = c("t-test", "wilcoxon", "t-test_overestim_var", "logreg"),
                                 corr_method = c("benjamini-hochberg", "bonferroni")) {
  
  library(reticulate)
  test <- match.arg(test)
  corr_method <- match.arg(corr_method)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  sc$tl$rank_genes_groups(scanpyObject, 
                          groupby = colDataName, 
                          groups = group1,
                          reference = group2,
                          method = test, 
                          n_genes = nGenes,
                          corr_method = corr_method,
                          layer = useAssay)
  
  py <- reticulate::py
  py$scanpyObject <- scanpyObject
  py_run_string("import pandas as pd")
  py_run_string(
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
  S4Vectors::metadata(inSCE)$scanpyMarkers <- markerGenesTable
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenes(sce, groups = '0')
#' }
#' @return plot object
#' @export
plotScanpyMarkerGenes <- function(inSCE,
                                  groups = NULL,
                                  nGenes = 10,
                                  nCols = 4,
                                  sharey = FALSE){
  
  if(is.null(inSCE@metadata[["findMarkerScanpyObject"]])){
    stop("marker genes not found. Kindly run runScanpyFindMarkers function first")
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesViolin(sce, groups = '0')
#' }
#' @return plot object
#' @export
plotScanpyMarkerGenesViolin <- function(inSCE,
                                        groups = NULL,
                                        features = NULL,
                                        nGenes = 10){
  
  if(is.null(inSCE@metadata["findMarkerScanpyObject"])){
    stop("marker genes not found. Kindly run runScanpyFindMarkers function first")
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesHeatmap(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
plotScanpyMarkerGenesHeatmap <- function(inSCE,
                                         groups = NULL,
                                         groupBy,
                                         nGenes = 10,
                                         features = NULL,
                                         log2fcThreshold = NULL){
  
  if(is.null(inSCE@metadata["findMarkerScanpyObject"])){
    stop("marker genes not found. Kindly run runScanpyFindMarkers function first")
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesDotPlot(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
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
  
  if(is.null(inSCE@metadata["findMarkerScanpyObject"])){
    stop("marker genes not found. Kindly run runScanpyFindMarkers function first")
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain_1" )
#' plotScanpyMarkerGenesMatrixPlot(sce, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
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
  
  if(is.null(inSCE@metadata["findMarkerScanpyObject"])){
    stop("marker genes not found. Kindly run runScanpyFindMarkers function first")
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyHeatmap(sce, features = markers, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
plotScanpyHeatmap <- function(inSCE,
                              features,
                              groupBy,
                              standardScale = 'var',
                              vmin = NULL,
                              vmax = NULL){
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyDotPlot(sce, features = markers, groupBy = 'Scanpy_louvain_1')
#' }
#' @return plot object
#' @export
plotScanpyDotPlot <- function(inSCE,
                              features,
                              groupBy,
                              standardScale = NULL,
                              title = '',
                              vmin = NULL,
                              vmax = NULL,
                              colorBarTitle = "Mean expression in group"){
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  
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


#' plotScanpyViolin
#' 
#' @param inSCE Input \code{SingleCellExperiment} object.
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
#' sce <- runScanpyPCA(sce, useAssay = "scanpyNormData")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' markers <- c("MALAT1" ,"RPS27" ,"CST3")
#' plotScanpyViolin(sce, features = markers, groupBy = "Scanpy_louvain_1")
#' }
#' @return plot object
#' @export
plotScanpyViolin <- function(inSCE,
                             features,
                             groupBy, 
                             xlabel = '', 
                             ylabel = NULL){
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  
  return(sc$pl$violin(scanpyObject,
                      keys = features,
                      groupby = groupBy,
                      xlabel = xlabel,
                      ylabel = ylabel))
  
}


