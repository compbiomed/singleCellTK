
# Helper/Wrapper Functions ---

#' .getComponentNames
#' Creates a list of PC/IC components to populate the picker for PC/IC heatmap
#' generation
#' @param maxComponents Number of components to return for the picker
#' @param component Which component to use. Choices are \code{PC} or \code{IC}.
#' @return List of component names (appended with \code{PC} or \code{IC})
.getComponentNames <- function(maxComponents, component = c("PC", "IC")){
  componentNames <- list()
  for (i in seq(maxComponents)) {
    componentNames[i] <- paste0(component, i)
  }
  return(componentNames)
}

#' .addSeuratToMetaDataSCE
#' Adds the input seurat object to the metadata slot of the input sce object
#' (after removing the data matrices)
#' @param inSCE (sce) object to which seurat object should be added in the
#' metadata slot (copy to)
#' @param seuratObject seurat object which should be added to the metadata slot
#' of sce object (copy from)
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' seurat object in its metadata slot (excluding data matrices)
.addSeuratToMetaDataSCE <- function(inSCE, seuratObject) {
  seuratObject@assays$RNA@counts <- methods::new("dgCMatrix")
  seuratObject@assays$RNA@data <- methods::new("dgCMatrix")
  seuratObject@assays$RNA@scale.data <- matrix()
  inSCE@metadata$seurat$obj <- seuratObject
  return(inSCE)
}

#' .computeSignificantPC
#' Computes the significant principal components from an input sce object (must
#' contain pca slot) using stdev
#' @param inSCE (sce) object with pca computed
#' @return A numerical value indicating how many number of components are
#' considered significant
.computeSignificantPC <- function(inSCE) {
  seuratObject <- convertSCEToSeurat(inSCE)
  max_components <- 0
  currentDistance <- 0
  previousDistance <- 0
  for (i in seq(seuratObject[["pca"]]@stdev)[-1]-1) {
    currentDistance <- abs(seuratObject[["pca"]]@stdev[i + 1] -
                             seuratObject[["pca"]]@stdev[i])
    if (abs(currentDistance - previousDistance) > 0.01) {
      previousDistance <- currentDistance
      max_components <- i
    }
    else{
      break
    }
  }
  return(max_components)
}

#' seuratNormalizeData
#' Wrapper for NormalizeData() function from seurat library
#' Normalizes the sce object according to the input parameters
#' @param inSCE (sce) object to normalize
#' @param useAssay Assay containing raw counts to use for normalization.
#' @param normAssayName Name of new assay containing normalized data. Default
#' \code{seuratNormData}.
#' @param normalizationMethod selected normalization method. Default
#' \code{"LogNormalize"}.
#' @param scaleFactor numeric value that represents the scaling factor. Default
#' \code{10000}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
seuratNormalizeData <- function(inSCE, useAssay,
                                normAssayName = "seuratNormData",
                                normalizationMethod = "LogNormalize",
                                scaleFactor = 10000, verbose = TRUE) {
  if(missing(useAssay)){
    useAssay <- SummarizedExperiment::assayNames(inSCE)[1]
    message("'useAssay' parameter missing. Using the first available assay ",
            "instead: '", useAssay, "'")
  }
  seuratObject <- Seurat::NormalizeData(convertSCEToSeurat(inSCE, useAssay),
                                        normalization.method = normalizationMethod,
                                        scale.factor = scaleFactor,
                                        verbose = verbose)
  inSCE <- .updateAssaySCE(inSCE, seuratObject, normAssayName, "data")
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  inSCE@metadata$seurat$normAssay <- normAssayName
  inSCE <- expSetDataTag(inSCE = inSCE, assayType = "normalized",
                         assays = normAssayName)
  return(inSCE)
}

#' seuratScaleData
#' Scales the input sce object according to the input parameters
#' @param inSCE (sce) object to scale
#' @param useAssay Assay containing normalized counts to scale.
#' @param scaledAssayName Name of new assay containing scaled data. Default
#' \code{seuratScaledData}.
#' @param model selected model to use for scaling data. Default \code{"linear"}.
#' @param scale boolean if data should be scaled or not. Default \code{TRUE}.
#' @param center boolean if data should be centered or not. Default \code{TRUE}
#' @param scaleMax maximum numeric value to return for scaled data. Default
#' \code{10}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' }
#' @return Scaled \code{SingleCellExperiment} object
#' @export
seuratScaleData <- function(inSCE, useAssay = "seuratNormData",
                            scaledAssayName = "seuratScaledData",
                            model = "linear", scale = TRUE, center = TRUE,
                            scaleMax = 10, verbose = TRUE) {
  seuratObject <- convertSCEToSeurat(inSCE, useAssay)
  seuratObject <- Seurat::ScaleData(seuratObject,
                                    features = rownames(seuratObject),
                                    model.use = model, do.scale = scale,
                                    do.center = center,
                                    scale.max = as.double(scaleMax),
                                    verbose = verbose)
  inSCE <- .updateAssaySCE(inSCE, seuratObject, scaledAssayName, "scale.data")
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  return(inSCE)
}

#' seuratFindHVG
#' Find highly variable genes and store in the input sce object
#' @param inSCE (sce) object to compute highly variable genes from and to store
#' back to it
#' @param useAssay Specify the name of the assay to use for computation
#'  of variable genes. It is recommended to use a raw counts assay with the 
#'  `vst` method and normalized assay with all other methods. Default
#'  is \code{"counts"}. 
#' @param hvgMethod selected method to use for computation of highly variable
#'  genes. One of 'vst', 'dispersion', or 'mean.var.plot'. Default method 
#'  is `vst` which uses the raw counts. All other methods use normalized counts.
#' @param hvgNumber numeric value of how many genes to select as highly
#' variable. Default \code{2000}
#' @param altExp Logical value indicating if the input object is an
#' altExperiment. Default \code{FALSE}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object with highly variable genes
#' computation stored
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
seuratFindHVG <- function(inSCE, useAssay = "counts",
                          hvgMethod = "vst", hvgNumber = 2000, altExp = FALSE,
                          verbose = TRUE) {
  
  if(hvgMethod == "vst"){
    seuratObject <- convertSCEToSeurat(inSCE, countsAssay = useAssay)
  }
  else{
    seuratObject <- convertSCEToSeurat(inSCE, normAssay = useAssay)
  }
  
  seuratObject <- Seurat::FindVariableFeatures(seuratObject,
                                               selection.method = hvgMethod,
                                               nfeatures = hvgNumber,
                                               verbose = verbose)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  if (hvgMethod == "vst") {
    if(!altExp){
      rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.variance.standardized
      rowData(inSCE)$seurat_variableFeatures_vst_mean <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.mean
    }
    else{
      # remove this part of code when updating to ExperimentSubset and add the
      # above code in if clause as the complete code
      altExpRows <- match(rownames(inSCE), rownames(methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features))
      rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.variance.standardized[altExpRows]
      rowData(inSCE)$seurat_variableFeatures_vst_mean <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.mean[altExpRows]
    }
  } else if (hvgMethod == "dispersion") {
    rowData(inSCE)$seurat_variableFeatures_dispersion_dispersion <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.dispersion
    rowData(inSCE)$seurat_variableFeatures_dispersion_dispersionScaled <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.dispersion.scaled
    rowData(inSCE)$seurat_variableFeatures_dispersion_mean <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.mean
  }
  else if (hvgMethod == "mean.var.plot") {
    rowData(inSCE)$seurat_variableFeatures_mvp_dispersion <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.dispersion
    rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.dispersion.scaled
    rowData(inSCE)$seurat_variableFeatures_mvp_mean <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$mvp.mean
  }
  return(inSCE)
}

#' seuratPCA
#' Computes PCA on the input sce object and stores the calculated principal
#' components within the sce object
#' @param inSCE (sce) object on which to compute PCA
#' @param useAssay Assay containing scaled counts to use in PCA.
#' @param reducedDimName Name of new reducedDims object containing Seurat PCA.
#' Default \code{seuratPCA}.
#' @param nPCs numeric value of how many components to compute. Default
#' \code{20}.
#' @param features Specify the feature names or rownames which should be used
#'  for computation of PCA. Default is \code{NULL} which will use the previously
#'  stored variable features.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' computed principal components
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratPCA <- function(inSCE, useAssay = "seuratScaledData",
                      reducedDimName = "seuratPCA", nPCs = 20, features = NULL,  verbose = TRUE) {
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  
  if(length(Seurat::VariableFeatures(seuratObject)) == 0
     && is.null(features)){
    features <- rownames(inSCE)
  }
  
  seuratObject <- Seurat::RunPCA(seuratObject,
                                 npcs = as.double(nPCs), verbose = verbose, features = features)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

  temp <- seuratObject@reductions$pca@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#' seuratICA
#' Computes ICA on the input sce object and stores the calculated independent
#' components within the sce object
#' @param inSCE (sce) object on which to compute ICA
#' @param useAssay Assay containing scaled counts to use in ICA.
#' @param reducedDimName Name of new reducedDims object containing Seurat ICA
#' Default \code{seuratICA}.
#' @param features Specify the feature names or rownames which should be used
#'  for computation of ICA. Default is \code{NULL} which will use the previously
#'  stored variable features.
#' @param nics Number of independent components to compute. Default \code{20}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratICA(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the
#' computed independent components
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratICA <- function(inSCE, useAssay,
                      reducedDimName = "seuratICA", features = NULL, nics = 20) {
  
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  
  if(length(Seurat::VariableFeatures(seuratObject)) == 0
     && is.null(features)){
    features <- rownames(inSCE)
  }
  
  seuratObject <- Seurat::RunICA(seuratObject,
                                 nics = as.double(nics), features = features, verbose = FALSE)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

  temp <- seuratObject@reductions$ica@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#' seuratComputeJackStraw
#' Compute jackstraw plot and store the computations in the input sce object
#' @param inSCE (sce) object on which to compute and store jackstraw plot
#' @param useAssay Assay containing scaled counts to use in JackStraw
#' calculation.
#' @param dims Number of components to test in Jackstraw. If \code{NULL}, then
#' all components are used. Default \code{NULL}.
#' @param numReplicate Numeric value indicating the number of replicate
#' samplings to perform.
#'  Default value is \code{100}.
#' @param propFreq Numeric value indicating the proportion of data to randomly
#' permute for each replicate.
#'  Default value is \code{0.025}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' sce <- seuratComputeJackStraw(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object with jackstraw
#' computations stored in it
#' @export
seuratComputeJackStraw <- function(inSCE, useAssay, dims = NULL,
                                   numReplicate = 100, propFreq = 0.025,
                                   externalReduction = NULL) {
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  if(!is.null(externalReduction)){
    #convert (_) to (-) as required by Seurat

    rownames(externalReduction@cell.embeddings) <- .convertToHyphen(rownames(externalReduction@cell.embeddings))
    seuratObject <- Seurat::FindVariableFeatures(seuratObject)
    seuratObject <- Seurat::ScaleData(seuratObject)
    seuratObject@reductions <- list(pca = externalReduction)
    seuratObject@reductions$pca@feature.loadings <- seuratObject@reductions$pca@feature.loadings[match(rownames(Seurat::GetAssayData(seuratObject, assay = "RNA", slot = "scale.data")), rownames(seuratObject@reductions$pca@feature.loadings)),]
    if(any(is.na(seuratObject@reductions$pca@feature.loadings))){
      seuratObject@reductions$pca@feature.loadings <- stats::na.omit(seuratObject@reductions$pca@feature.loadings)
    }
    seuratObject@commands$RunPCA.RNA <- seuratObject@commands$ScaleData.RNA
    seuratObject@commands$RunPCA.RNA@params$rev.pca <- FALSE
    seuratObject@commands$RunPCA.RNA@params$weight.by.var <- TRUE
    }
  if(is.null(seuratObject@reductions[["pca"]])) {
    stop("'seuratPCA' must be run before JackStraw can be computed.")
  }
  if(is.null(dims)) {
    dims <- ncol(seuratObject@reductions[["pca"]])
  }
  seuratObject <- Seurat::JackStraw(seuratObject, dims = as.double(dims),
                                    num.replicate = numReplicate,
                                    prop.freq = propFreq)
  seuratObject <- Seurat::ScoreJackStraw(seuratObject, dims = seq(dims))
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  return(inSCE)
}

#' seuratJackStrawPlot
#' Computes the plot object for jackstraw plot from the pca slot in the input
#' sce object
#' @param inSCE (sce) object from which to compute the jackstraw plot (pca
#' should be computed)
#' @param dims Number of components to plot in Jackstraw. If \code{NULL}, then
#' all components are plotted Default \code{NULL}.
#' @param xmax X-axis maximum on each QQ plot. Default \code{0.1}.
#' @param ymax Y-axis maximum on each QQ plot. Default \code{0.3}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' sce <- seuratComputeJackStraw(sce, useAssay = "counts")
#' seuratJackStrawPlot(sce)
#' }
#' @return plot object
#' @export
seuratJackStrawPlot <- function(inSCE, dims = NULL, xmax = 0.1, ymax = 0.3,
                                externalReduction = NULL) {
  seuratObject <- convertSCEToSeurat(inSCE)
  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
  }
  if(is.null(seuratObject@reductions[["pca"]])) {
    stop("'seuratPCA' must be run before JackStraw can be computed.")
  }
  if(is.null(dims)) {
    dims <- ncol(seuratObject@reductions[["pca"]])
  }
  return(Seurat::JackStrawPlot(seuratObject, dims = seq(dims),
                               xmax = xmax, ymax = ymax))
}

#' seuratPlotHVG
#' Plot highly variable genes from input sce object (must have highly variable
#' genes computations stored)
#' @param inSCE (sce) object that contains the highly variable genes
#' computations
#' @param labelPoints Numeric value indicating the number of top genes that
#' should be labeled.
#'  Default is \code{0}, which will not label any point.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' seuratPlotHVG(sce)
#' }
#' @return plot object
#' @export
seuratPlotHVG <- function(inSCE, labelPoints = 0) {
  seuratObject <- convertSCEToSeurat(inSCE)
  plot <- Seurat::VariableFeaturePlot(seuratObject)
  plot$labels$colour <- "Variable"
  if(requireNamespace("stringr", quietly = TRUE)){
    plot$data$colors <- stringr::str_to_title(plot$data$colors)
  }
  if(labelPoints > 0){
    plot <- Seurat::LabelPoints(plot,
                        points = Seurat::VariableFeatures(
                          object = seuratObject)[seq(labelPoints)])
  }
  return(plot)
}

#' seuratReductionPlot
#' Plots the selected dimensionality reduction method
#' @param inSCE (sce) object which has the selected dimensionality reduction
#' algorithm already computed and stored
#' @param useReduction Dimentionality reduction to plot. One of "pca", "ica",
#' "tsne", or "umap". Default \code{"umap"}.
#' @param showLegend Select if legends and labels should be shown on the output
#' plot or not. Either "TRUE" or "FALSE". Default \code{FALSE}.
#' @param groupBy Specify a colData column name that be used for grouping.
#' Default is \code{NULL}.
#' @param splitBy Specify a colData column name that be used for splitting the
#' output plot. Default is \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' seuratReductionPlot(sce, useReductionPlot = "pca")}
#' @return plot object
#' @export
seuratReductionPlot <- function(inSCE, useReduction = c("pca", "ica",
                                                        "tsne", "umap"),
                                showLegend = FALSE, groupBy = NULL,
                                splitBy = NULL) {
  seuratObject <- convertSCEToSeurat(inSCE)
  if(!is.null(seuratObject@meta.data$seurat_cluster)){
    seuratObject@meta.data <- seuratObject@meta.data[, "seurat_clusters", drop=FALSE]
  }
  else{
    seuratObject@meta.data <- data.frame()
  }

  if(showLegend){
    if(!is.null(seuratObject@meta.data$seurat_clusters)){
      Seurat::Idents(seuratObject) <- seuratObject@meta.data$seurat_clusters
      seuratObject@meta.data <- data.frame()
    }
  }

  if(!is.null(groupBy)){
    seuratObject[[groupBy]] <- colData(inSCE)[[groupBy]]
  }

  if(!is.null(splitBy)){
    seuratObject[[splitBy]] <- colData(inSCE)[[splitBy]]
  }

  if(showLegend){
    plot <- Seurat::DimPlot(
      object = seuratObject,
      reduction = useReduction,
      group.by = groupBy,
      split.by = splitBy,
      label = TRUE)
  }
  else{
    plot <- Seurat::DimPlot(
      object = seuratObject,
      reduction = useReduction,
      group.by = groupBy,
      split.by = splitBy,
      label = FALSE) + Seurat::NoLegend()
  }

  if ("ident" %in% names(plot$data) &&
      "seurat_clusters" %in% names(seuratObject@meta.data)) {
    plot$data$ident <- seuratObject@meta.data$seurat_clusters
  }
  return(plot)
}


#' seuratFindClusters
#' Computes the clusters from the input sce object and stores them back in sce
#' object
#' @param inSCE (sce) object from which clusters should be computed and stored
#' in
#' @param useAssay Assay containing scaled counts to use for clustering.
#' @param useReduction Reduction method to use for computing clusters. One of
#' "pca" or "ica". Default \code{"pca"}.
#' @param dims numeric value of how many components to use for computing
#' clusters. Default \code{10}.
#' @param algorithm selected algorithm to compute clusters. One of "louvain",
#' "multilevel", or "SLM". Use \code{louvain} for "original Louvain algorithm"
#' and \code{multilevel} for "Louvain algorithm with multilevel refinement".
#' Default \code{louvain}.
#' @param groupSingletons boolean if singletons should be grouped together or
#' not. Default \code{TRUE}.
#' @param resolution Set the resolution parameter to find larger (value above 1)
#' or smaller (value below 1) number of communities. Default \code{0.8}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' sce <- seuratFindClusters(sce, useAssay = "counts")
#' }
#' @return Updated sce object which now contains the computed clusters
#' @export
seuratFindClusters <- function(
  inSCE,
  useAssay = "seuratScaledData",
  useReduction = c("pca", "ica"),
  dims = 10,
  algorithm = c("louvain", "multilevel", "SLM"),
  groupSingletons = TRUE,
  resolution = 0.8,
  externalReduction = NULL,
  verbose = TRUE) {

  algorithm <- match.arg(algorithm)
  useReduction <- match.arg(useReduction)

  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)

  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
    useReduction <- "pca"
  }

  seuratObject <- Seurat::FindNeighbors(seuratObject, reduction = useReduction,
                                        dims = seq(dims), verbose = verbose)
  no_algorithm <- 1
  if (algorithm == "louvain") {
    no_algorithm = 1
  } else if (algorithm == "multilevel") {
    no_algorithm = 2
  } else if (algorithm == "SLM") {
    no_algorithm = 3
  }
  tempSeuratObject <- seuratObject
  tempSeuratObject@meta.data <- data.frame()
  tempSeuratObject <- Seurat::FindClusters(tempSeuratObject,
                                           algorithm = no_algorithm,
                                           group.singletons = groupSingletons,
                                           resolution = resolution,
                                           verbose = verbose)
  seuratObject@meta.data$seurat_clusters <- tempSeuratObject@meta.data$seurat_clusters
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  colData(inSCE)[[paste0("Seurat","_",algorithm,"_","Resolution",resolution)]] <- seuratObject@meta.data$seurat_clusters
  S4Vectors::metadata(inSCE)$seurat$clusterName <- paste0("Seurat", "_",
                                                          algorithm, "_",
                                                          "Resolution",
                                                          resolution)
  return(inSCE)
}

#' seuratRunTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back
#' into the sce object
#' @param inSCE (sce) object on which to compute the tSNE
#' @param useReduction selected reduction algorithm to use for computing tSNE.
#' One of "pca" or "ica". Default \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Seurat tSNE
#' Default \code{seuratTSNE}.
#' @param dims Number of reduction components to use for tSNE computation.
#' Default \code{10}.
#' @param perplexity Adjust the perplexity tuneable parameter for the underlying
#' tSNE call. Default \code{30}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @return Updated sce object with tSNE computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratRunTSNE <- function(inSCE, useReduction = c("pca", "ica"),
                          reducedDimName = "seuratTSNE", dims = 10,
                          perplexity = 30, externalReduction = NULL) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE)
  
  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
    useReduction <- "pca"
  }
  
  seuratObject <- Seurat::RunTSNE(seuratObject, reduction = useReduction,
                                  dims = seq(dims), perplexity = perplexity,
                                  check_duplicates = FALSE)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  temp <- seuratObject@reductions$tsne@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#' seuratRunUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back
#' into the sce object
#' @param inSCE (sce) object on which to compute the UMAP
#' @param useReduction Reduction to use for computing UMAP. One of "pca" or
#' "ica". Default is \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Seurat UMAP
#' Default \code{seuratUMAP}.
#' @param dims Numerical value of how many reduction components to use for UMAP
#' computation. Default \code{10}.
#' @param minDist Sets the \code{"min.dist"} parameter to the underlying UMAP
#' call. See \link[Seurat]{RunUMAP} for more information. Default \code{0.3}.
#' @param nNeighbors Sets the \code{"n.neighbors"} parameter to the underlying
#' UMAP call. See \link[Seurat]{RunUMAP} for more information. Default
#' \code{30L}.
#' @param spread Sets the \code{"spread"} parameter to the underlying UMAP call.
#' See \link[Seurat]{RunUMAP} for more information. Default \code{1}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' sce <- seuratFindClusters(sce, useAssay = "counts")
#' sce <- seuratRunUMAP(sce, useReduction = "pca")
#' }
#' @return Updated sce object with UMAP computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratRunUMAP <- function(inSCE, useReduction = c("pca", "ica"),
                          reducedDimName = "seuratUMAP", dims = 10,
                          minDist = 0.3, nNeighbors = 30L, spread = 1,
                          externalReduction = NULL,
                          verbose = TRUE) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE)
  
  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
    useReduction <- "pca"
  }
  
  seuratObject <- Seurat::RunUMAP(seuratObject,
                                  reduction = useReduction,
                                  dims = seq(dims),
                                  min.dist = minDist,
                                  n.neighbors = nNeighbors,
                                  spread = spread,
                                  verbose = verbose)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

  temp <- seuratObject@reductions$umap@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#' .seuratGetVariableFeatures
#' Retrieves the requested number of variable feature names
#' @param inSCE (sce) object from which to extract the variable feature names
#' @param numberOfFeatures numerical value indicating how many feature names
#' should be retrieved (default is 100)
#' @return list() of variable feature names
.seuratGetVariableFeatures <- function(inSCE, numberOfFeatures) {
  seuratObject <- convertSCEToSeurat(inSCE)
  if (length(seuratObject@assays$RNA@var.features) > 0) {
    return(seuratObject@assays$RNA@var.features[seq(numberOfFeatures)])
  }
}

#' seuratElbowPlot
#' Computes the plot object for elbow plot from the pca slot in the input sce
#' object
#' @param inSCE (sce) object from which to compute the elbow plot (pca should
#' be computed)
#' @param significantPC Number of significant principal components to plot.
#' This is used to alter the color of the points for the corresponding PCs.
#' If \code{NULL}, all points will be the same color. Default \code{NULL}.
#' @param reduction Reduction to use for elbow plot generation. Either
#' \code{"pca"} or \code{"ica"}. Default \code{"pca"}.
#' @param ndims Number of components to use. Default \code{20}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @param interactive Logical value indicating if the returned object should
#'  be an interactive plotly object if \code{TRUE} or a ggplot object if
#'  set to \code{FALSE}. Default is \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' seuratElbowPlot(sce)
#' }
#' @return plot object
#' @export
seuratElbowPlot <- function(inSCE,
                            significantPC = NULL,
                            reduction = "pca",
                            ndims = 20,
                            externalReduction = NULL,
                            interactive = TRUE) {
  seuratObject <- convertSCEToSeurat(inSCE)
  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
    reduction <- "pca"
  }
  plot <- Seurat::ElbowPlot(seuratObject, reduction = reduction, ndims = ndims)
  if(!is.null(significantPC)){
    if(significantPC > ndims){
      significantPC <- ndims
    }
    plot$data$Significant <- c(rep("Yes", significantPC),
                               rep("No", length(rownames(plot$data)) - significantPC))
    plot <- ggplot2::ggplot(data = plot$data,
                            ggplot2::aes(x = plot$data$dims,
                                         y = plot$data$stdev,
                                         color = plot$data$Significant)) +
      ggplot2::geom_point()
  }
  plot$labels$x <- "PC"
  plot$labels$y <- "Standard Deviation"
  plot$labels$colour <- "Significant"

  if(interactive){
    hoverText <- paste("Dimension:", plot$data$dims, "\nStandard Deviation:",
                       round(plot$data$stdev, 1), "\nIs Significant?",
                       plot$data$Significant)
    significant <- plot$data$Significant
    if(length(unique(significant))>1){
      plot <- plotly::style(plot, text = hoverText[seq(which(significant == "No")[1])-1])
      plot <- plotly::style(plot, text = hoverText[which(significant == "No")[1]:length(significant)], traces = 1)
    }
    else{
      plot <- plotly::style(plot, text = hoverText)
    }
  }

  return(plot)
}

#' seuratComputeHeatmap
#' Computes the heatmap plot object from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute heatmap (pca should be
#' computed)
#' @param useAssay Assay containing scaled counts to use in heatmap.
#' @param useReduction Reduction method to use for computing clusters. One of
#' "pca" or "ica". Default \code{"pca"}.
#' @param dims Number of components to generate heatmap plot objects. If
#' \code{NULL}, a heatmap will be generated for all components. Default
#' \code{NULL}.
#' @param nfeatures Number of features to include in the heatmap. Default
#' \code{30}.
#' @param cells Numeric value indicating the number of top cells to plot.
#'  Default is \code{NULL} which indicates all cells.
#' @param ncol Numeric value indicating the number of columns to use for plot.
#'  Default is \code{NULL} which will automatically compute accordingly.
#' @param balanced Plot equal number of genes with positive and negative scores.
#'  Default is \code{TRUE}.
#' @param fast See \link[Seurat]{DimHeatmap} for more information. Default
#' \code{TRUE}.
#' @param combine See \link[Seurat]{DimHeatmap} for more information. Default
#' \code{TRUE}.
#' @param raster See \link[Seurat]{DimHeatmap} for more information. Default
#' \code{TRUE}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through
#' other libraries. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' heatmap <- seuratComputeHeatmap(sce, useAssay = "counts")
#' seuratHeatmapPlot(heatmap)
#' }
#' @return plot object
#' @export
seuratComputeHeatmap <- function(inSCE,
                                 useAssay,
                                 useReduction = c("pca", "ica"),
                                 dims = NULL,
                                 nfeatures = 30,
                                 cells = NULL,
                                 ncol = NULL,
                                 balanced = TRUE,
                                 fast = TRUE,
                                 combine = TRUE,
                                 raster = TRUE,
                                 externalReduction = NULL) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  if(!is.null(externalReduction)){
    seuratObject@reductions <- list(pca = externalReduction)
    useReduction <- "pca"
  }
  if(is.null(dims)) {
    dims <- ncol(seuratObject@reductions[[useReduction]])
  }
  return(Seurat::DimHeatmap(seuratObject,
                            dims = seq(dims),
                            nfeatures = nfeatures,
                            cells = cells,
                            reduction = useReduction,
                            ncol = ncol,
                            fast = fast,
                            combine = combine,
                            raster = raster,
                            balanced = balanced))
}

#' seuratHeatmapPlot
#' Modifies the heatmap plot object so it contains specified number of heatmaps
#' in a single plot
#' @param plotObject plot object computed from seuratComputeHeatmap() function
#' @param dims numerical value of how many heatmaps to draw (default is 0)
#' @param ncol numerical value indicating that in how many columns should the
#' heatmaps be distrbuted (default is 2)
#' @param labels list() of labels to draw on heatmaps
#' @return modified plot object
#' @export
seuratHeatmapPlot <- function(plotObject, dims, ncol, labels) {
  componentsToPlot <- as.integer(gsub("[^0-9.]", "", labels))
  return(cowplot::plot_grid(plotlist = plotObject[c(componentsToPlot)],
                            ncol = ncol, labels = labels))
}

#' .updateAssaySCE
#' Update/Modify/Add an assay in the provided SingleCellExperiment object from
#' a Seurat object
#' @param inSCE Input SingleCellExperiment object
#' @param seuratObject Input Seurat object
#' @param assaySlotSCE Selected assay to update in the input
#' SingleCellExperiment object
#' @param seuratDataSlot Selected data slot from the Seurat object. Default
#' \code{"counts"}.
#' @param seuratAssaySlot Selected assay from Seurat object. Default
#' \code{"RNA"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  data from Seurat object appended to the \link{assay} slot.
#' @importFrom SummarizedExperiment assay<-
.updateAssaySCE <- function(inSCE, seuratObject, assaySlotSCE,
                            seuratDataSlot = "counts",
                            seuratAssaySlot = "RNA") {
  assay(inSCE, assaySlotSCE) <- NULL
  temp.matrix <- methods::slot(Seurat::GetAssay(seuratObject, seuratAssaySlot),
                               seuratDataSlot)
  rownames(temp.matrix) <- rownames(inSCE)
  colnames(temp.matrix) <- colnames(inSCE)
  assay(inSCE, assaySlotSCE) <- temp.matrix
  return(inSCE)
}

#' convertSeuratToSCE
#' Converts the input seurat object to a sce object
#' @param seuratObject Input Seurat object
#' @param normAssayName Name of assay to store the normalized data. Default
#' \code{"seuratNormData"}.
#' @param scaledAssayName Name of assay to store the scaled data. Default
#' \code{"seuratScaledData"}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' seurat <- convertSCEToSeurat(sce)
#' sce <- convertSeuratToSCE(seurat)
#' @return \code{SingleCellExperiment} output object
#' @export
convertSeuratToSCE <- function(seuratObject, normAssayName = "seuratNormData",
                               scaledAssayName = "seuratScaledData") {
  inSCE <- SingleCellExperiment(
    assays = list(counts = seuratObject@assays[[1]]@counts),
    colData = seuratObject@meta.data)
  
  assay(inSCE, normAssayName) <- methods::slot(seuratObject@assays$RNA, "data")
  if (length(methods::slot(seuratObject, "assays")[["RNA"]]@scale.data) > 0) {
    assay(inSCE, scaledAssayName) <- methods::slot(seuratObject@assays$RNA,
                                                   "scale.data")
  }
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  return(inSCE)
}

#' convertSCEToSeurat
#' Converts sce object to seurat while retaining all assays and metadata
#' @param inSCE A \code{SingleCellExperiment} object to convert to a Seurat
#'  object.
#' @param countsAssay Which assay to use from sce object for raw counts.
#'  Default \code{NULL}.
#' @param normAssay Which assay to use from sce object for normalized data.
#'  Default \code{NULL}.
#' @param scaledAssay Which assay to use from sce object for scaled data.
#'  Default \code{NULL}.
#' @param copyColData Boolean. Whether copy 'colData' of SCE object to
#'  the 'meta.data' of Seurat object. Default \code{FALSE}.
#' @param copyReducedDim Boolean. Whether copy 'reducedDims' of the SCE
#'  object to the 'reductions' of Seurat object. Default \code{FALSE}.
#' @param copyDecontX Boolean. Whether copy 'decontXcounts' assay of the
#'  SCE object to the 'assays' of Seurat object. Default \code{TRUE}.
#' @param pcaReducedDim Specify a character value indicating the name of
#'  the reducedDim to store as default pca computation in the output seurat
#'  object. Default is \code{NULL} which will not store any reducedDim as the
#'  default pca. This will only work when \code{copyReducedDim} parameter is
#'  set to \code{TRUE}.
#' @param icaReducedDim Specify a character value indicating the name of
#'  the reducedDim to store as default ica computation in the output seurat
#'  object. Default is \code{NULL} which will not store any reducedDim as the
#'  default ica. This will only work when \code{copyReducedDim} parameter is
#'  set to \code{TRUE}.
#' @param tsneReducedDim Specify a character value indicating the name of
#'  the reducedDim to store as default tsne computation in the output seurat
#'  object. Default is \code{NULL} which will not store any reducedDim as the
#'  default tsne. This will only work when \code{copyReducedDim} parameter is
#'  set to \code{TRUE}.
#' @param umapReducedDim Specify a character value indicating the name of
#'  the reducedDim to store as default umap computation in the output seurat
#'  object. Default is \code{NULL} which will not store any reducedDim as the
#'  default umap. This will only work when \code{copyReducedDim} parameter is
#'  set to \code{TRUE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' seurat <- convertSCEToSeurat(sce)
#' @return Updated seurat object that contains all data from the input sce
#' object
#' @export
#' @importFrom SummarizedExperiment assay assays
convertSCEToSeurat <- function(inSCE, countsAssay = NULL, normAssay = NULL,
                               scaledAssay = NULL, copyColData = FALSE,
                               copyReducedDim = FALSE, copyDecontX = FALSE,
                               pcaReducedDim = NULL, icaReducedDim = NULL,
                               tsneReducedDim = NULL, umapReducedDim = NULL) {

  .checkSCEValidity(inSCE)

  if(!is.null(countsAssay) && !(countsAssay %in% expDataNames(inSCE))) {
    stop(paste0("'", countsAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
  }
  if(!is.null(normAssay) && !(normAssay %in% expDataNames(inSCE))) {
    stop(paste0("'", normAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
  }
  if(!is.null(scaledAssay) && !(scaledAssay %in% expDataNames(inSCE))) {
    stop(paste0("'", scaledAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
  }
  if(!is.null(pcaReducedDim) && !(pcaReducedDim %in% reducedDimNames(inSCE))){
    stop(paste0("'", pcaReducedDim, "' not found in the list of reducedDims: ",
                paste(reducedDimNames(inSCE), collapse=",")))
  }
  if(!is.null(icaReducedDim) && !(icaReducedDim %in% reducedDimNames(inSCE))){
    stop(paste0("'", icaReducedDim, "' not found in the list of reducedDims: ",
                paste(reducedDimNames(inSCE), collapse=",")))
  }
  if(!is.null(tsneReducedDim) && !(tsneReducedDim %in% reducedDimNames(inSCE))){
    stop(paste0("'", tsneReducedDim, "' not found in the list of reducedDims: ",
                paste(reducedDimNames(inSCE), collapse=",")))
  }
  if(!is.null(umapReducedDim) && !(umapReducedDim %in% reducedDimNames(inSCE))){
    stop(paste0("'", umapReducedDim, "' not found in the list of reducedDims: ",
                paste(reducedDimNames(inSCE), collapse=",")))
  }

  # Seurat has a particular way of modifying row/colnames
  # Save row/colnames in metadata
  seuratRowNames <- gsub("_", "-", rownames(inSCE))
  seuratColNames <- gsub("_", "-", colnames(inSCE))
  inSCE@metadata$seurat$colNames <- seuratColNames
  inSCE@metadata$seurat$rowNames <- seuratRowNames

  # Create Seurat object and Set counts assay
  # If no counts assay is supplied, the first assay is used
  if (!is.null(countsAssay) && countsAssay %in% names(assays(inSCE))) {
    temp <- .convertToMatrix(assay(inSCE, countsAssay))
  } else {
    temp <- .convertToMatrix(assays(inSCE)[[1]])
  }
  rownames(temp) <- seuratRowNames
  colnames(temp) <- seuratColNames
  seuratObject <- Seurat::CreateSeuratObject(counts = temp)

  # Set normalized assay
  if (!is.null(normAssay) && normAssay %in% names(assays(inSCE))) {
    tempMatrix <- .convertToMatrix(assay(inSCE, normAssay))
    if(inherits(tempMatrix, "dgeMatrix")){
      tempMatrix <- methods::as(tempMatrix, "dgCMatrix")
    }
    seuratObject@assays$RNA@data <- tempMatrix
    rownames(seuratObject@assays$RNA@data) <- seuratRowNames
    colnames(seuratObject@assays$RNA@data) <- seuratColNames
  }

  # Set Scaled Assay
  if (!is.null(scaledAssay) && scaledAssay %in% names(assays(inSCE))) {
    seuratObject@assays$RNA@scale.data <- as.matrix(assay(inSCE, scaledAssay))
    rownames(seuratObject@assays$RNA@scale.data) <- seuratRowNames
    colnames(seuratObject@assays$RNA@scale.data) <- seuratColNames
  }

  if (!is.null(inSCE@metadata$seurat$obj)) {
    if(length(inSCE@metadata$seurat$obj@assays$RNA@var.features) > 0) {
      seuratObject@assays$RNA@var.features <- inSCE@metadata$seurat$obj@assays$RNA@var.features
    }
    if (!is.null(inSCE@metadata$seurat$obj@reductions$pca)) {
      seuratObject@reductions$pca <- inSCE@metadata$seurat$obj@reductions$pca
    }
    if (!is.null(inSCE@metadata$seurat$obj@assays$RNA@meta.features)) {
      seuratObject@assays$RNA@meta.features <- inSCE@metadata$seurat$obj@assays$RNA@meta.features
    }
    if (!is.null(inSCE@metadata$seurat$obj@reductions$ica)) {
      seuratObject@reductions$ica <- inSCE@metadata$seurat$obj@reductions$ica
    }
    if (!is.null(inSCE@metadata$seurat$obj@reductions$tsne)) {
      seuratObject@reductions$tsne <- inSCE@metadata$seurat$obj@reductions$tsne
    }
    if (!is.null(inSCE@metadata$seurat$obj@reductions$umap)) {
      seuratObject@reductions$umap <- inSCE@metadata$seurat$obj@reductions$umap
    }
    if (!is.null(inSCE@metadata$seurat$obj@meta.data)) {
      seuratObject@meta.data <- inSCE@metadata$seurat$obj@meta.data
    }
    if (!is.null(inSCE@metadata$seurat$obj@commands)) {
      seuratObject@commands <- inSCE@metadata$seurat$obj@commands
    }
  }

  # Set colData from inSCE object if required
  if (!is.null(colData(inSCE)) && copyColData) {
    seuratObject@meta.data <- cbind(seuratObject@meta.data, colData(inSCE))
  }

  # Set additional reducedDims from inSCE object if required
  if (length(SingleCellExperiment::reducedDims(inSCE)) > 0 && copyReducedDim) {
    for (redc in SingleCellExperiment::reducedDimNames(inSCE)) {
      reDim <- SingleCellExperiment::reducedDim(inSCE, redc)
      colnames(reDim) <- paste0(redc, "_", seq_len(length(colnames(reDim))))
      rownames(reDim) <- gsub('_', '-', rownames(reDim))
      key <-  gsub('_', '', redc)
      seuratObject@reductions[[redc]] <- Seurat::CreateDimReducObject(embeddings = reDim,
                                                                key = paste0(key, "_"), assay = "RNA")
    }

    availReducedDims <- c("pca", "ica", "tsne", "umap")
    for(i in availReducedDims){
      if(!is.null(eval(parse(text = paste0(i, "ReducedDim"))))){
        seuratObject@reductions[[i]] <- seuratObject@reductions[[eval(parse(text = paste0(i, "ReducedDim")))]]
        seuratObject@reductions[[eval(parse(text = paste0(i, "ReducedDim")))]] <- NULL
        message("'", eval(parse(text = paste0(i, "ReducedDim"))), "' reducedDim from input SCE object saved to the default ", i, " slot of seurat object.")
      }
    }
  }

  # Set 'decontXCounts' assay to seurat object if required
  if ("decontXcounts" %in% SummarizedExperiment::assayNames(inSCE) && copyDecontX) {
    decontM <- SummarizedExperiment::assay(inSCE, "decontXcounts")
    colnames(decontM) <- colnames(seuratObject)
    rownames(decontM) <- gsub('_', '-', rownames(decontM))
    seuratObject[["decontXcounts"]] <- Seurat::CreateAssayObject(counts = .convertToMatrix(decontM))
  }
  
  # Ensuring that colnames from input SCE converted to Seurat object are same in the Seurat metadata slot
  rownames(seuratObject@meta.data) <- colnames(seuratObject)

  return(seuratObject)
}

#' seuratSCTransform
#' Runs the \link[Seurat]{SCTransform} function to transform/normalize the input
#' data
#' @param inSCE Input SingleCellExperiment object
#' @param normAssayName Name for the output data assay. Default
#' \code{"SCTCounts"}.
#' @param useAssay Name for the input data assay. Default \code{"counts"}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @return Updated SingleCellExperiment object containing the transformed data
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- seuratSCTransform(mouseBrainSubsetSCE)
seuratSCTransform <- function(inSCE, normAssayName = "SCTCounts",
                              useAssay = "counts", verbose = TRUE) {
  seuratObject <- base::suppressWarnings(Seurat::SCTransform(
    object = convertSCEToSeurat(inSCE, useAssay),
    assay = "RNA",
    new.assay.name = "SCTransform",
    do.correct.umi = FALSE,
    verbose = verbose))
  inSCE <- .updateAssaySCE(inSCE = inSCE, seuratObject = seuratObject,
                           assaySlotSCE = normAssayName,
                           seuratDataSlot = "data",
                           seuratAssaySlot = "SCTransform")
  inSCE <- expSetDataTag(inSCE = inSCE, assayType = "normalized",
                         assays = normAssayName)
  return(inSCE)
}


#' .seuratInvalidate
#' Removes seurat data from the input SingleCellExperiment object specified by
#' the task in the Seurat workflow.
#' @param inSCE Input \code{SingleCellExperiment} object to remove Seurat data
#' from.
#' @param scaleData Remove scaled data from seurat. Default \code{TRUE}.
#' @param varFeatures Remove variable features from seurat. Default \code{TRUE}.
#' @param PCA Remove PCA from seurat. Default \code{TRUE}.
#' @param ICA Remove ICA from seurat. Default \code{TRUE}.
#' @param tSNE Remove tSNE from seurat. Default \code{TRUE}.
#' @param UMAP Remove UMAP from seurat. Default \code{TRUE}.
#' @param clusters Remove clusters from seurat. Default \code{TRUE}.
#' @return Updated SingleCellExperiment object containing the Seurat object in
#' the metadata slot with the data removed
#' @importFrom SummarizedExperiment assay<-
.seuratInvalidate <- function(inSCE, scaleData = TRUE, varFeatures = TRUE,
                              PCA = TRUE, ICA = TRUE, tSNE = TRUE, UMAP = TRUE,
                              clusters = TRUE){
  if(scaleData){
    assay(inSCE, "seuratScaledData") <- NULL
  }
  if(varFeatures){
    methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@var.features <- logical()
    methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features <- data.frame(row.names = make.unique(gsub("_", "-", rownames(inSCE))))
    inSCE@metadata$seurat$heatmap_pca <- NULL
  }
  if(PCA){
    inSCE@metadata$seurat$obj@reductions$pca <- NULL
  }
  if(ICA){
    inSCE@metadata$seurat$obj@reductions$ica <- NULL
  }
  if(tSNE){
    inSCE@metadata$seurat$obj@reductions$tsne <- NULL
  }
  if(UMAP){
    inSCE@metadata$seurat$obj@reductions$umap <- NULL
  }
  if(clusters){
    inSCE@metadata$seurat$obj@meta.data$seurat_clusters <- NULL
  }
  return(inSCE)
}


#' seuratIntegration
#' A wrapper function to Seurat Batch-Correction/Integration workflow.
#' @param inSCE Input \code{SingleCellExperiment} object that contains the assay
#' to batch-correct.
#' @param useAssay Assay to batch-correct.
#' @param batch Batch variable from \code{colData} slot of
#' \code{SingleCellExperiment} object.
#' @param newAssayName Assay name for the batch-corrected output assay.
#' @param kAnchor Number of neighbours to use for finding the anchors in the
#' \link[Seurat]{FindIntegrationAnchors} function.
#' @param kFilter Number of neighbours to use for filtering the anchors in the
#' \link[Seurat]{FindIntegrationAnchors} function.
#' @param kWeight Number of neighbours to use when weigthing the anchors in the
#' \link[Seurat]{IntegrateData} function.
#' @param ndims Number of dimensions to use. Default \code{10}.
#'
#' @return A \code{SingleCellExperiment} object that contains the
#' batch-corrected assay inside the \code{altExp} slot of the object
#' @export
seuratIntegration <- function(inSCE, useAssay = "counts", batch,
                              newAssayName = "SeuratIntegratedAssay", kAnchor,
                              kFilter, kWeight, ndims = 10){
  if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
    stop(paste(useAssay, "not found in the input object assays"))
  }
  if(is.null(batch)){
    stop("batch variable must be provided for batch-correction")
  }
  if(kAnchor == 0 || kFilter == 0 || kWeight == 0){
    stop("kAnchor, kFilter or kWeight cannot be zero. Please input correct ",
         "parameters.")
  }

  #create seurat object
  seuratObject <- convertSCEToSeurat(inSCE, useAssay)
  rownames(seuratObject@meta.data) <- gsub("_", "-",
                                           rownames(seuratObject@meta.data))
  seuratObject@meta.data <- cbind(seuratObject@meta.data, colData(inSCE))

  #split seurat object by batch variable
  seurat.list <- Seurat::SplitObject(seuratObject, split.by = batch)
  seurat.list <- seurat.list[c(unique(seuratObject@meta.data[[batch]]))]

  #find anchors
  seurat.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat.list,
                                                   dims = seq(ndims),
                                                   k.anchor = kAnchor,
                                                   k.filter = kFilter)
  seurat.integrated <- Seurat::IntegrateData(anchorset = seurat.anchors,
                                             dims = seq(ndims),
                                             k.weight = kWeight)
  #store results back in altExp slot of sce object
  altExp(inSCE, newAssayName) <- SingleCellExperiment(list(counts = Seurat::GetAssayData(seurat.integrated@assays$integrated, "data")))
  SummarizedExperiment::assayNames(altExp(inSCE,newAssayName)) <- newAssayName
  # remove this if counts in above line set to altExp

  # store back colData from sce into the altExp slot
  colData(altExp(inSCE, newAssayName))<- colData(inSCE)

  #counts <- assay(altExp(inSCE, newAssayName), "altExp")
  # remove NA values from counts and replace with zero so can be used properly
  # by dgCMatrix
  counts <- assay(altExp(inSCE, newAssayName), newAssayName)
  counts[is.na(counts)] <- 0

  #store back counts
  #assay(altExp(inSCE, newAssayName), "altExp") <- counts
  assay(altExp(inSCE, newAssayName), newAssayName) <- counts

  return(inSCE)
}

#' seuratFindMarkers
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param cells1 A \code{list} of sample names included in group1.
#' @param cells2 A \code{list} of sample names included in group2.
#' @param group1 Name of group1.
#' @param group2 Name of group2.
#' @param allGroup Name of all groups.
#' @param conserved Logical value indicating if markers conserved between two
#' groups should be identified. Default is \code{FALSE}.
#' @param test Test to use for DE. Default \code{"wilcox"}.
#' @param onlyPos Logical value indicating if only positive markers should be
#' returned.
#' @param minPCT Numeric value indicating the minimum fraction of min.pct
#'  cells in which genes are detected. Default is \code{0.1}.
#' @param threshUse Numeric value indicating the logFC threshold value on
#'  which on average, at least X-fold difference (log-scale) between the
#'  two groups of cells exists. Default is \code{0.25}.
#' @param verbose Logical value indicating if informative messages should
#'  be displayed. Default is \code{TRUE}.
#' @return A \code{SingleCellExperiment} object that contains marker genes
#' populated in a data.frame stored inside metadata slot.
#' @export
seuratFindMarkers <- function(
  inSCE, cells1 = NULL, cells2 = NULL, group1 = NULL, group2 = NULL,
  allGroup = NULL, conserved = FALSE, test = "wilcox", onlyPos = FALSE,
  minPCT = 0.1, threshUse = 0.25, verbose = TRUE){
  seuratObject <- convertSCEToSeurat(inSCE)
  markerGenes <- NULL
  if(is.null(allGroup)
     && (!is.null(group1) && !is.null(group2))){
    #convert (_) to (-) as required by Seurat
    cells1 <- .convertToHyphen(cells1)
    cells2 <- .convertToHyphen(cells2)
    Seurat::Idents(seuratObject, cells = cells1) <- group1
    Seurat::Idents(seuratObject, cells = cells2) <- group2
    markerGenes <- NULL
    if(!conserved){
      markerGenes <- Seurat::FindMarkers(
        object = seuratObject,
        ident.1 = group1,
        ident.2 = group2,
        test.use = test,
        only.pos = onlyPos
        )
    }
    else{
      seuratObject[["groups"]] <- Seurat::Idents(seuratObject)
      markerGenes <- .findConservedMarkers(
        object = seuratObject,
        ident.1 = group1,
        ident.2 = group2,
        grouping.var = "groups",
        test.use = test,
        only.pos = onlyPos,
        cells1 = cells1,
        cells2 = cells2)
    }
    markerGenes$cluster1 <- group1
    markerGenes$cluster2 <- group2
    gene.id <- rownames(markerGenes)
    markerGenes <- cbind(gene.id, markerGenes)
  }
  else if(!is.null(allGroup)
           && (is.null(group1) && is.null(group2))){
    Seurat::Idents(seuratObject,
                   cells = colnames(seuratObject)) <- colData(inSCE)[[allGroup]]
    markerGenes <- Seurat::FindAllMarkers(
      seuratObject,
      test.use = test,
      only.pos = onlyPos,
      logfc.threshold = threshUse,
      min.pct = minPCT,
      verbose = verbose)
    gene.id <- markerGenes$gene
    markerGenes <- cbind(gene.id, markerGenes)
    markerGenes$gene <- NULL
    # grp <- unique(colData(inSCE)[[allGroup]])
    # clust <- as.integer(unique(Seurat::Idents(seuratObject)))
    # for(i in seq(length(clust))){
    #   levels(markerGenes$cluster)[clust[i]] <- grp[i]
    # }
    colnames(markerGenes)[which(colnames(markerGenes) == "cluster")] <- "cluster1"
    markerGenes$cluster2 <- rep("all", nrow(markerGenes))
  }
  else if(is.null(allGroup)
          && (is.null(group1) && is.null(group2))){
    Seurat::Idents(
      seuratObject,
      cells = colnames(seuratObject)) <- colData(inSCE)[[S4Vectors::metadata(inSCE)$seurat$clusterName]]
    markerGenes <- Seurat::FindAllMarkers(
      seuratObject,
      test.use = test,
      only.pos = onlyPos,
      logfc.threshold = threshUse,
      min.pct = minPCT,
      verbose = verbose)
    gene.id <- markerGenes$gene
    markerGenes <- cbind(gene.id, markerGenes)
    markerGenes$gene <- NULL
    # grp <- unique(colData(inSCE)[[S4Vectors::metadata(inSCE)$seurat$clusterName]])
    # clust <- as.integer(unique(Seurat::Idents(seuratObject)))
    # for(i in seq(length(clust))){
    #   levels(markerGenes$cluster)[clust[i]] <- grp[i]
    # }
  }
  rownames(markerGenes) <- NULL
  S4Vectors::metadata(inSCE)$seuratMarkers <- markerGenes
  return(inSCE)
}

#' Compute and plot visualizations for marker genes
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param scaledAssayName Specify the name of the scaled assay stored in the
#' input object.
#' @param plotType Specify the type of the plot to compute. Options are limited
#' to "ridge", "violing", "feature", "dot" and "heatmap".
#' @param features Specify the features to compute the plot against.
#' @param groupVariable Specify the column name from the colData slot that
#' should be used as grouping variable.
#' @param splitBy Specify the column name from the colData slot that should be
#' used to split samples.
#'  Default is \code{NULL}.
#' @param cols Specify two colors to form a gradient between. Default is
#' \code{c("lightgrey", "blue")}.
#' @param ncol Visualizations will be adjusted in "ncol" number of columns.
#'  Default is \code{1}.
#'
#' @return Plot object
#' @export
seuratGenePlot <- function(inSCE,
                           scaledAssayName = "seuratScaledData",
                           plotType,
                           features,
                           groupVariable,
                           splitBy = NULL,
                           cols = c("lightgrey", "blue"),
                           ncol = 1){
  #setup seurat object and the corresponding groups
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = scaledAssayName)
  indices <- list()
  cells <- list()
  groups <- unique(colData(inSCE)[[groupVariable]])
  for(i in seq(length(groups))){
    indices[[i]] <- which(colData(inSCE)[[groupVariable]] == groups[i],
                          arr.ind = TRUE)
    cells[[i]] <- colnames(inSCE)[indices[[i]]]
    cells[[i]] <- .convertToHyphen(cells[[i]])
    Seurat::Idents(seuratObject, cells = cells[[i]]) <- groups[i]
  }

  if(!is.null(splitBy)){
    seuratObject[[splitBy]] <- colData(inSCE)[[splitBy]]
  }

  #plot required visualization
  if(plotType == "ridge"){
    return(Seurat::RidgePlot(
      seuratObject,
      features = features,
      ncol = ncol))
  }
  else if(plotType == "violin"){
    return(Seurat::VlnPlot(
      seuratObject,
      features = features,
      ncol = ncol,
      split.by = splitBy))
  }
  else if(plotType == "feature"){
    return(Seurat::FeaturePlot(
      seuratObject,
      features = features,
      cols = cols,
      ncol = ncol,
      split.by = splitBy))
  }
  else if(plotType == "dot"){
    return(Seurat::DotPlot(
      seuratObject,
      features = features,
      split.by = splitBy))
  }
  else if(plotType == "heatmap"){
    return(Seurat::DoHeatmap(seuratObject, features = features))
  }
}

.findConservedMarkers <- function(object,
                                  ident.1,
                                  ident.2,
                                  grouping.var = "groups",
                                  test.use,
                                  only.pos,
                                  cells1,
                                  cells2){
  meta.method <- metap::minimump
  verbose <- TRUE
  slot <- "data"
  assay <- "RNA"
  marker.test <- list()

  object.var <- Seurat::FetchData(object = object, vars = grouping.var)
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)

  marker.test[[1]] <- Seurat::FindMarkers(
    object = object,
    assay = assay,
    slot = slot,
    ident.1 = levels.split[1],
    ident.2 = levels.split[2],
  )

  marker.test[[2]] <- Seurat::FindMarkers(
    object = object,
    assay = assay,
    slot = slot,
    ident.1 = levels.split[2],
    ident.2 = levels.split[3],
  )

  marker.test[[3]] <- Seurat::FindMarkers(
    object = object,
    assay = assay,
    slot = slot,
    ident.1 = levels.split[3],
    ident.2 = levels.split[1],
  )

  names(x = marker.test)[1] <- levels.split[1]
  names(x = marker.test)[2] <- levels.split[2]
  names(x = marker.test)[3] <- levels.split[3]

  marker.test <- Filter(f = Negate(f = is.null), x = marker.test)
  genes.conserved <- Reduce(
    f = intersect,
    x = lapply(
      X = marker.test,
      FUN = function(x) {
        return(rownames(x = x))
      }
    )
  )
  markers.conserved <- list()
  for (i in seq_len(length(x = marker.test))) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, ]
    colnames(x = markers.conserved[[i]]) <- paste(
      names(x = marker.test)[i],
      colnames(x = markers.conserved[[i]]),
      sep = "_"
    )
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  pval.codes <- colnames(x = markers.combined)[grepl(pattern = "*_p_val$", x = colnames(x = markers.combined))]
  if (length(x = pval.codes) > 1) {
    markers.combined$max_pval <- apply(
      X = markers.combined[, pval.codes, drop = FALSE],
      MARGIN = 1,
      FUN = max
    )
    combined.pval <- data.frame(cp = apply(
      X = markers.combined[, pval.codes, drop = FALSE],
      MARGIN = 1,
      FUN = function(x) {
        return(meta.method(x)$p)
      }
    ))
    meta.method.name <- as.character(x = formals()$meta.method)
    if (length(x = meta.method.name) == 3) {
      meta.method.name <- meta.method.name[3]
    }
    colnames(x = combined.pval) <- paste0(meta.method.name, "_p_val")
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined <- markers.combined[order(markers.combined[, paste0(meta.method.name, "_p_val")]), ]
  }
  lfcCol <- colnames(markers.combined)[grep(paste0(ident.1, "_avg"), colnames(markers.combined))]
  markers.combined <- markers.combined[, c(
    paste0(ident.1, "_p_val"),
    lfcCol,
    paste0(ident.1, "_pct.1"),
    paste0(ident.1, "_pct.2"),
    paste0("_p_val"))]
  colnames(markers.combined) <- gsub(pattern = paste0(ident.1, "_"),
                                     replacement = "",
                                     x = colnames(markers.combined))
  colnames(markers.combined) <- c(colnames(markers.combined)[-length(colnames(markers.combined))], "p_val_adj")
  return(markers.combined)
}

#' Get variable feature names after running seuratFindHVG function
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#'
#' @return A list of variable feature names.
#' @export
seuratVariableFeatures <- function(inSCE){
  if(!is.null(S4Vectors::metadata(inSCE)$seurat$obj)){
    return(Seurat::VariableFeatures(S4Vectors::metadata(inSCE)$seurat$obj))
  }
  else{
    return(NULL)
  }
}
