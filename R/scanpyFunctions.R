
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
#' @param normAssayName Name of new assay containing normalized data. Default
#' \code{scanpyNormData}.
#' @param normalizationMethod selected normalization method. Default
#' \code{"LogNormalize"}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
runScanpyNormalizeData <- function(inSCE,
                                   useAssay,
                                   normAssayName = "scanpyNormData",
                                   normalizationMethod = "LogNormalize") {
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
  if(normalizationMethod == "LogNormalize"){
    # Total-count normalize (library-size correct) to 10,000 reads/cell
    sc$pp$normalize_per_cell(scanpyObject, counts_per_cell_after = 1e4)
    # log transform the data.
    sc$pp$log1p(scanpyObject)
  }
  else{
    scanpyObject <- sc$pp$normalize_total(scanpyObject, 
                                          target_sum=1e4, 
                                          inplace = FALSE)
  }
  
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
#' variable. Default \code{2000}
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
                             hvgNumber = 2000) {
  
  method <- match.arg(method)
  scanpyObject <- convertSCEToScanpy(inSCE, useAssay)
  sc$pp$highly_variable_genes(scanpyObject, flavor = method)
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
  } else if (method == "cell_ranger") {
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
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "counts")
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
#' @param algorithm selected method to use for computation of pca. One of \code{'arpack'}, \code{'randomized'}, \code{'auto'} or \code{'lobpcg'}.
#' Default \code{"arpack"}.
#' @param use_highly_variable boolean value of whether to use highly variable genes only. By default uses them if they have been determined beforehand. 
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runScanpyNormalizeData(sce, useAssay = "counts")
#' sce <- runScanpyFindHVG(sce, useAssay = "counts")
#' sce <- runScanpyScaleData(sce, useAssay = "counts")
#' sce <- runScanpyPCA(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' computed principal components
#' @export
#' @importFrom SingleCellExperiment reducedDim<- rowSubset
#' @importFrom S4Vectors metadata<-
#' 
runScanpyPCA <- function(inSCE,
                         useAssay = "scanpyScaledData",
                         reducedDimName = "scanpyPCA",
                         nPCs = 50L,
                         algorithm = c("arpack", "randomized", "auto", "lobpcg"),
                         use_highly_variable = TRUE){
  
  
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
#' @param useReduction Reduction to use for computing UMAP. Default is \code{"pca"}.
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
#' #' @param alpha Sets the \code{"alpha"} parameter to the underlying UMAP call.
#' Default \code{1}.
#' #' @param gamma Sets the \code{"gamma"} parameter to the underlying UMAP call.
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
                          dims = 10L,  #ncomponents
                          minDist = 0.5,
                          nNeighbors = 15L, #n.neighbors
                          spread = 1,
                          alpha=1.0, 
                          gamma=1.0, 
                          externalReduction = NULL) {
  
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
                          dims = 10L,#ncomponents 
                          perplexity = 15L,
                          externalReduction = NULL){
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE, X_name = useAssay)
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
#' "leiden", or "dendogram". Default \code{louvain}.
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
runSeuratFindClusters <- function(inSCE,
                                  useAssay = "scanpyScaledData",
                                  useReduction = "scanpyPCA",
                                  nNeighbors = 15L,
                                  dims = 2L,
                                  algorithm = c("louvain", "leiden", "dendrogram"),
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
  } else if (algorithm == "dendrogram") {
    sc$tl$dendrogram(adata = scanpyObject,
                     n_pcs = dims,
                     key_added = colDataName,
                     use_rep = useReduction,
                     cor_method = cor_method,
                     inplace = TRUE)
  }
  
  
  inSCE@metadata$scanpy$obj <- scanpyObject
  colData(inSCE)[[paste0("Scanpy", "_", algorithm)]] <-
    as.factor(unlist(scanpyObject$obs[algorithm]))
  S4Vectors::metadata(inSCE)$scanpy[algorithm] <- paste0("Scanpy", "_",
                                                         algorithm)
  return(inSCE)
}


runScanpyFindMarkers <- function(inSCE,
                                 nGenes,
                                 clusterName,
                                 group1 = "all",
                                 group2 = "rest",
                                 test = c("logreg", "t-test", "wilcoxon", "t-test_overestim_var"),
                                 corr_method = c("benjamini-hochberg", "bonferroni")) {
  
  test <- match.arg(test)
  corr_method <- match.arg(corr_method)
  
  scanpyObject <- zellkonverter::SCE2AnnData(sce = inSCE)
  markerGenes <- NULL
  
  sc$tl$rank_genes_groups(scanpyObject, 
                          groupby = clusterName, 
                          groups = group1,
                          reference = group2,
                          method = test, 
                          n_genes = nGenes,
                          method = test,
                          corr_method = corr_method)
  ##after this i want to access scanpyObject$uns['rank_genes_groups]['names'] 
  #to store these marker genes names within sce object 
  rownames(markerGenes) <- NULL
  S4Vectors::metadata(inSCE)$scanpyMarkers <- markerGenes
  return(inSCE)
}

################################################
###EDITING ALREADY EXISTING FUNCTION WITHIN SCTK
################################################
.dfFromHVGMetric <- function(inSCE,
                             method = c("vst", "mean.var.plot", "dispersion",
                                        "modelGeneVar", "seurat_v3", 
                                        "cell_ranger", "seurat")) {
  method <- match.arg(method)
  df <- data.frame(featureNames = rownames(inSCE))
  if (method == "vst") {
    m <- "seurat_variableFeatures_vst_mean"
    v_rank <- "seurat_variableFeatures_vst_varianceStandardized"
    v_plot <- "seurat_variableFeatures_vst_varianceStandardized"
    if (is.null(rowData(inSCE)[[v_rank]]) || is.null(rowData(inSCE)[[m]])) {
      stop("Seurat vst metric not found in inSCE. ",
           "Run `runSeuratFindHVG()` with 'vst' method before ",
           "using this function!")
    }
  } else if (method == "dispersion") {
    m <- "seurat_variableFeatures_dispersion_mean"
    v_rank <- "seurat_variableFeatures_dispersion_dispersion"
    v_plot <- "seurat_variableFeatures_dispersion_dispersionScaled"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Seurat dispersion metric not found in inSCE. ",
           "Run `runSeuratFindHVG()` with 'dispersion' method ",
           "before using this function!")
    }
  } else if (method == "seurat_v3") {
    m <- "scanpy_variableFeatures_seuratv3_mean"
    v_rank <- "scanpy_variableFeatures_seuratv3_variances"
    v_plot <- "scanpy_variableFeatures_seuratv3_variancesScaled"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Scanpy variance metric not found in inSCE. ",
           "Run `runScanpyFindHVG()` with 'seurat_v3' method ",
           "before using this function!")
    }
  }else if (method == "cell_ranger") {
    m <- "scanpy_variableFeatures_cr_mean"
    v_rank <- "scanpy_variableFeatures_cr_dispersion"
    v_plot <- "scanpy_variableFeatures_cr_dispersionScaled"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Scanpy dispersion metric not found in inSCE. ",
           "Run `runScanpyFindHVG()` with 'cell_ranger' method ",
           "before using this function!")
    }
  }else if (method == "seurat") {
    m <- "scanpy_variableFeatures_seurat_mean"
    v_rank <- "scanpy_variableFeatures_seurat_dispersion"
    v_plot <- "scanpy_variableFeatures_seurat_dispersionScaled"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Scanpy dispersion metric not found in inSCE. ",
           "Run `runScanpyFindHVG()` with 'dispersion' method ",
           "before using this function!")
    }
  }else if (method == "modelGeneVar") {
    m <- "scran_modelGeneVar_mean"
    v_rank <- "scran_modelGeneVar_bio"
    v_plot <- "scran_modelGeneVar_totalVariance"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Scran modelGeneVar metric not found in inSCE. Run ",
           "`runModelGeneVar()` before using this function!")
    }
  } else if (method == "mean.var.plot") {
    m <- "seurat_variableFeatures_mvp_mean"
    v_rank <- "seurat_variableFeatures_mvp_dispersion"
    v_plot <- "seurat_variableFeatures_mvp_dispersionScaled"
    if (is.null(rowData(inSCE)[[v_rank]]) ||
        is.null(rowData(inSCE)[[m]]) ||
        is.null(rowData(inSCE)[[v_plot]])) {
      stop("Seurat mean.var.plot metric not found in inSCE. ",
           "Run `runSeuratFindHVG()` with ",
           "'mean.var.plot' method before using this function!")
    }
  }
  df$mean <- rowData(inSCE)[[m]]
  df$v_rank <- rowData(inSCE)[[v_rank]]
  df$v_plot <- rowData(inSCE)[[v_plot]]
  return(df)
}

getTopHVG <- function(inSCE,
                      method = c("vst", "dispersion",
                                 "mean.var.plot", "modelGeneVar", "seurat", 
                                 "seurat_v3", "cell_ranger"),
                      hvgNumber = 2000,
                      useFeatureSubset = NULL,
                      featureDisplay = metadata(inSCE)$featureDisplay) {
  method <- match.arg(method)
  topGenes <- character()
  if (!is.null(useFeatureSubset)) {
    topGenes <- .parseUseFeatureSubset(inSCE, useFeatureSubset,
                                       altExpObj = NULL, returnType = "cha")
  } else {
    metrics <- .dfFromHVGMetric(inSCE, method)
    metrics <- metrics[order(-metrics$v_rank),]
    metrics <- metrics[metrics$v_rank > 0, ]
    if (method == "mean.var.plot") {
      means.use <- (metrics$mean > 0.1) & (metrics$mean < 8)
      dispersions.use <- (metrics$v_plot > 1) & (metrics$v_plot < Inf)
      metrics <- metrics[means.use & dispersions.use,]
    }
    hvgNumber <- min(hvgNumber, nrow(metrics))
    topGenes <- as.character(metrics$featureNames)[seq_len(hvgNumber)]
  }
  topGenes <- stats::na.omit(topGenes)
  if (!is.null(featureDisplay)) {
    geneIdx <- featureIndex(topGenes, inSCE)
    topGenes <- rowData(inSCE)[[featureDisplay]][geneIdx]
  }
  return(topGenes)
}
##############################
##############################