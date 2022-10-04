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

#' runScanpyNormalizeData
#' Wrapper for NormalizeData() function from scanpy library
#' Normalizes the sce object according to the input parameters
#' @param inSCE (sce) object to normalize
#' @param useAssay Assay containing raw counts to use for normalization.
#' @param countsPerCellAfter If None, after normalization, each cell has a total
#' count equal to the median of the counts_per_cell before normalization.
#' @param countsPerCell Precomputed counts per cell.
#' @param minCount Cells with counts less than min_counts are filtered out 
#' during normalization.
#' @param normAssayName Name of new assay containing normalized data. Default
#' \code{scanpyNormData}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
runScanpyNormalizeData <- function(inSCE,
                                   useAssay,
                                   countsPerCellAfter = 1e4,
                                   countsPerCell = NULL,
                                   minCount = 0,
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
  # Total-count normalize (library-size correct) to 10,000 reads/cell
  sc$pp$normalize_per_cell(scanpyObject, 
                           counts_per_cell_after = countsPerCellAfter, 
                           counts_per_cell = countsPerCell,
                           min_counts = minCount)
  # log transform the data.
  sc$pp$log1p(scanpyObject)
  
  inSCE <-
    .updateAssaySCEFromScanpy(inSCE, scanpyObject, normAssayName)
  inSCE@metadata$scanpy$obj <- scanpyObject
  inSCE@metadata$scanpy$normAssay <- normAssayName
  inSCE <- expSetDataTag(inSCE = inSCE,
                         assayType = "normalized",
                         assays = normAssayName)
  return(inSCE)
}

#' runScanpyFindHVG
#' Find highly variable genes and store in the input sce object
#' @param inSCE (sce) object to compute highly variable genes from and to store
#' back to it
#' @param useAssay Specify the name of the assay to use for computation
#'  of variable genes. It is recommended to use a raw counts assay with the
#' @param method selected method to use for computation of highly variable
#' genes. One of \code{'seurat'}, \code{'cell_ranger'}, or \code{'seurat_v3'}.
#' Default \code{"seurat"}.
#' @param hvgNumber numeric value of how many genes to select as highly
#' variable. Default \code{NULL}
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
#' sce <- runScanpyFindHVG(sce)
#' @return Updated \code{SingleCellExperiment} object with highly variable genes
#' computation stored
#' \code{\link{getTopHVG}}, \code{\link{plotTopHVG}}
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors metadata<-
runScanpyFindHVG <- function(inSCE,
                             useAssay = "counts",
                             method = c("seurat", "cell_ranger", "seurat_v3"),
                             hvgNumber = NULL,
                             minMean = 0.0125,
                             maxMean = 3,
                             minDisp = 0.5,
                             maxDisp = Inf) {
  
  method <- match.arg(method)
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
  sc$pp$scale(scanpyObject, max_value=10)
  sc$pp$highly_variable_genes(scanpyObject, 
                              flavor = method,
                              n_top_genes = hvgNumber,
                              min_mean = minMean,
                              max_mean = maxMean,
                              min_disp = minDisp,
                              max_disp = maxDisp)
  #options
  #sc["pp"]["highly_variable_genes"](scanpyObject, min_mean, max_mean, min_disp)
  #inSCE <- convertScanpytoSCE(scanpyObject)
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  if (method == "seurat") {
    rowData(inSCE)$scanpy_variableFeatures_seurat_dispersion <-
      inSCE@metadata$scanpy$obj["var"][["dispersions"]]
    rowData(inSCE)$scanpy_variableFeatures_seurat_dispersionScaled <-
      inSCE@metadata$scanpy$obj["var"][["dispersions_norm"]]    
    rowData(inSCE)$scanpy_variableFeatures_seurat_mean <-
      inSCE@metadata$scanpy$obj["var"][["means"]]  
    
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
  #for this approach it is required that filering of cells must be done 
  #initially so that duplicates can be dropped
  else if (method == "cell_ranger") {
    rowData(inSCE)$scanpy_variableFeatures_cr_dispersion <-
      inSCE@metadata$scanpy$obj["var"][["dispersions"]]
    rowData(inSCE)$scanpy_variableFeatures_cr_dispersionScaled <-
      inSCE@metadata$scanpy$obj["var"][["dispersions_norm"]]    
    rowData(inSCE)$scanpy_variableFeatures_cr_mean <-
      inSCE@metadata$scanpy$obj["var"][["means"]]  
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
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_variances <-
      inSCE@metadata$scanpy$obj["var"][["variances"]]
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_variancesScaled <-
      inSCE@metadata$scanpy$obj["var"][["variances_norm"]]    
    rowData(inSCE)$scanpy_variableFeatures_seuratv3_mean <-
      inSCE@metadata$scanpy$obj["var"][["means"]]  
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
  inSCE@metadata$scanpy$obj <- scanpyObject
  inSCE@metadata$scanpy$scaledAssay <- scaledAssayName
  inSCE <- expSetDataTag(inSCE = inSCE,
                         assayType = "scaled",
                         assays = scaledAssayName)
  
  return(inSCE)
}


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
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
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
  
  sc$tl$pca(scanpyObject, svd_solver= algorithm, n_comps = nPCs)
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  
  temp <- scanpyObject$obsm['X_pca']
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  
  return(inSCE)
}


#' runScanpyUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back
#' into the sce object
#' @param inSCE (sce) object on which to compute the UMAP
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
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts")
#' sce <- runScanpyUMAP(sce, useReduction = "scanpyPCA")
#' }
#' @return Updated sce object with UMAP computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
runScanpyUMAP <- function(inSCE,
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
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReduction <- "pca"
  } 
  sc$pp$neighbors(scanpyObject, 
                  n_neighbors = nNeighbors, 
                  n_pcs = dims,
                  use_rep = useReduction)
  sc$tl$umap(scanpyObject, 
             n_components = dims,
             min_dist = minDist,
             alpha = alpha, 
             gamma = gamma,
             spread = spread)
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  
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
                          useReduction = "scanpyPCA", 
                          reducedDimName = "scanpyTSNE",
                          dims = 10L,
                          perplexity = 15L,
                          externalReduction = NULL){
  
  params <- as.list(environment())
  params$inSCE <- NULL
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  if (!is.null(externalReduction)) {
    scanpyObject$obsm <- list(pca = externalReduction)
    useReduction <- "pca"
  } 
  
  sc$tl$tsne(scanpyObject, 
             n_pcs = dims,
             use_rep = useReduction, 
             perplexity = perplexity)
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  
  temp <- scanpyObject$obsm['X_tsne']
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp
  metadata(inSCE)$sctk$runDimReduce$reddim[[reducedDimName]] <- params
  
  return(inSCE)
}



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
#' @param resolution A parameter value controlling the coarseness of the 
#' clustering. Higher values lead to more clusters Default \code{1}.
#' @param colDataName colName to store the result in annData object
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
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "counts")
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
                                  resolution,
                                  colDataName,
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
  
  if(missing(colDataName)){
    colDataName = algorithm
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
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  colData(inSCE)[[paste0("Scanpy", "_", algorithm)]] <-
    as.factor(unlist(scanpyObject$obs[algorithm]))
  S4Vectors::metadata(inSCE)$scanpy[algorithm] <- paste0("Scanpy", "_",
                                                         algorithm)
  return(inSCE)
}


#' runScanpyFindMarkers
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param nGenes The number of genes that appear in the returned tables. 
#' Defaults to all genes.
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
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' sce <- runScanpyFindClusters(sce, useAssay = "counts", algorithm = "louvain")
#' sce <- runScanpyFindMarkers(sce, colDataName = "Scanpy_louvain" )
#' }
#' @return A \code{SingleCellExperiment} object that contains marker genes
#' populated in a data.frame stored inside metadata slot.
#' @export
runScanpyFindMarkers <- function(inSCE,
                                 nGenes = NULL,
                                 colDataName,
                                 group1 = "all",
                                 group2 = "rest",
                                 test = c("t-test", "wilcoxon", "t-test_overestim_var", "logreg"),
                                 corr_method = c("benjamini-hochberg", "bonferroni")) {
  
  test <- match.arg(test)
  corr_method <- match.arg(corr_method)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  sc$pp$normalize_per_cell(scanpyObject, 
                           counts_per_cell_after = 1e4, 
                           counts_per_cell = NULL,
                           min_counts = 0)
  # log transform the data.
  sc$pp$log1p(scanpyObject)
  sc$tl$rank_genes_groups(scanpyObject, 
                          groupby = colDataName, 
                          groups = group1,
                          reference = group2,
                          method = test, 
                          n_genes = nGenes,
                          corr_method = corr_method)
  
  py_run_string("import pandas as pd")
  py_run_string(
    "names = pd.DataFrame(r.scanpyObject.uns['rank_genes_groups']['names'])", 
    convert = TRUE)
  py_run_string(
    "logFoldChanges = 
    pd.DataFrame(r.scanpyObject.uns['rank_genes_groups']['logfoldchanges'])", 
    convert = TRUE)
  py_run_string(
    "pvals_adj = 
    pd.DataFrame(r.scanpyObject.uns['rank_genes_groups']['pvals_adj'])", 
    convert = TRUE)
  py_run_string(
    "scores = 
    pd.DataFrame(r.scanpyObject.uns['rank_genes_groups']['scores'])",
    convert = TRUE)
  
  markerGenesNames <- utils::stack(py$names)
  colnames(markerGenesNames) <- c("Gene","findMarker_cluster")
  Log2_FC <- unlist(py$logFoldChanges)
  Pvalue <- unlist(py$pvals_adj)
  zscore <- unlist(py$scores)
  
  
  markerGenesTable <- data.frame()
  markerGenesTable <- cbind(markerGenesNames, Log2_FC, Pvalue, zscore)
  
  S4Vectors::metadata(inSCE)$scanpyMarkers <- markerGenesTable
  return(inSCE)
}
