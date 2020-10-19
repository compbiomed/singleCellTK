
# Helper/Wrapper Functions ---

#' .getComponentNames
#' Creates a list of PC/IC components to populate the picker for PC/IC heatmap generation
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
#' Adds the input seurat object to the metadata slot of the input sce object (after removing the data matrices)
#' @param inSCE (sce) object to which seurat object should be added in the metadata slot (copy to)
#' @param seuratObject seurat object which should be added to the metadata slot of sce object (copy from)
#' @return Updated \code{SingleCellExperiment} object which now contains the seurat object in its metadata slot (excluding data matrices)
.addSeuratToMetaDataSCE <- function(inSCE, seuratObject) {
    seuratObject@assays$RNA@counts <- methods::new("dgCMatrix")
    seuratObject@assays$RNA@data <- methods::new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    inSCE@metadata$seurat$obj <- seuratObject
    return(inSCE)
}

#' .computeSignificantPC
#' Computes the significant principal components from an input sce object (must contain pca slot) using stdev
#' @param inSCE (sce) object with pca computed
#' @return A numerical value indicating how many number of components are considered significant
.computeSignificantPC <- function(inSCE) {
    seuratObject <- convertSCEToSeurat(inSCE)
    max_components <- 0
    for (i in seq(seuratObject[["pca"]]@stdev)[-1]-1) {
        if (abs(seuratObject[["pca"]]@stdev[i + 1] - seuratObject[["pca"]]@stdev[i]) > 0.1) {
            max_components <- i
        }
    }
    return(max_components)
}

#' seuratNormalizeData
#' Wrapper for NormalizeData() function from seurat library
#' Normalizes the sce object according to the input parameters
#' @param inSCE (sce) object to normalize
#' @param useAssay Assay containing raw counts to use for normalization.
#' @param normAssayName Name of new assay containing normalized data. Default \code{seuratNormData}.
#' @param normalizationMethod selected normalization method. Default \code{"LogNormalize"}.
#' @param scaleFactor numeric value that represents the scaling factor. Default \code{10000}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' }
#' @return Normalized \code{SingleCellExperiment} object
#' @export
seuratNormalizeData <- function(inSCE, useAssay, normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000) {
    seuratObject <- Seurat::NormalizeData(convertSCEToSeurat(inSCE, useAssay), normalization.method = normalizationMethod, scale.factor = scaleFactor, verbose = FALSE)
    inSCE <- .updateAssaySCE(inSCE, seuratObject, normAssayName, "data")
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    inSCE@metadata$seurat$normAssay <- normAssayName
    return(inSCE)
}

#' seuratScaleData
#' Scales the input sce object according to the input parameters
#' @param inSCE (sce) object to scale
#' @param useAssay Assay containing normalized counts to scale.
#' @param scaledAssayName Name of new assay containing scaled data. Default \code{seuratScaledData}.
#' @param model selected model to use for scaling data. Default \code{"linear"}.
#' @param scale boolean if data should be scaled or not. Default \code{TRUE}.
#' @param center boolean if data should be centered or not. Default \code{TRUE}
#' @param scaleMax maximum numeric value to return for scaled data. Default \code{10}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' }
#' @return Scaled \code{SingleCellExperiment} object
#' @export
seuratScaleData <- function(inSCE, useAssay, scaledAssayName = "seuratScaledData", model = "linear", scale = TRUE, center = TRUE, scaleMax = 10) {
  seuratObject <- convertSCEToSeurat(inSCE, useAssay)
  seuratObject <- Seurat::ScaleData(seuratObject, features = rownames(seuratObject), model.use = model, do.scale = scale, do.center = center, scale.max = as.double(scaleMax), verbose = FALSE)
  inSCE <- .updateAssaySCE(inSCE, seuratObject, scaledAssayName, "scale.data")
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  return(inSCE)
}

#' seuratFindHVG
#' Find highly variable genes and store in the input sce object
#' @param inSCE (sce) object to compute highly variable genes from and to store back to it
#' @param useAssay Normalized assay inside the SCE object to use for hvg computation.
#' @param hvgMethod selected method to use for computation of highly variable genes. One of 'vst', 'dispersion', or 'mean.var.plot'. Default \code{"vst"}.
#' @param hvgNumber numeric value of how many genes to select as highly variable. Default \code{2000}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object with highly variable genes computation stored
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
seuratFindHVG <- function(inSCE, useAssay,  hvgMethod = "vst", hvgNumber = 2000) {
  seuratObject <- convertSCEToSeurat(inSCE, normAssay = useAssay)
  seuratObject <- Seurat::FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber, verbose = FALSE)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  if (hvgMethod == "vst") {
    rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.variance.standardized
    rowData(inSCE)$seurat_variableFeatures_vst_mean <- methods::slot(inSCE@metadata$seurat$obj, "assays")[["RNA"]]@meta.features$vst.mean
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
#' Computes PCA on the input sce object and stores the calculated principal components within the sce object
#' @param inSCE (sce) object on which to compute PCA
#' @param useAssay Assay containing scaled counts to use in PCA.
#' @param reducedDimName Name of new reducedDims object containing Seurat PCA. Default \code{seuratPCA}.
#' @param nPCs numeric value of how many components to compute. Default \code{20}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the computed principal components
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratPCA <- function(inSCE, useAssay, reducedDimName = "seuratPCA", nPCs = 20) {
    seuratObject <- Seurat::RunPCA(convertSCEToSeurat(inSCE, scaledAssay = useAssay), npcs = as.double(nPCs), verbose = FALSE)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

    temp <- seuratObject@reductions$pca@cell.embeddings
    rownames(temp) <- colnames(inSCE)
    reducedDim(inSCE, reducedDimName) <- temp

    return(inSCE)
}

#' seuratICA
#' Computes ICA on the input sce object and stores the calculated independent components within the sce object
#' @param inSCE (sce) object on which to compute ICA
#' @param useAssay Assay containing scaled counts to use in ICA.
#' @param reducedDimName Name of new reducedDims object containing Seurat ICA Default \code{seuratICA}.
#' @param nics Number of independent components to compute. Default \code{20}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratICA(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object which now contains the computed independent components
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratICA <- function(inSCE, useAssay, reducedDimName = "seuratICA", nics = 20) {
    seuratObject <- Seurat::RunICA(convertSCEToSeurat(inSCE, scaledAssay = useAssay), nics = as.double(nics), verbose = FALSE)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

    temp <- seuratObject@reductions$ica@cell.embeddings
    rownames(temp) <- colnames(inSCE)
    reducedDim(inSCE, reducedDimName) <- temp

    return(inSCE)
}

#' seuratComputeJackStraw
#' Compute jackstraw plot and store the computations in the input sce object
#' @param inSCE (sce) object on which to compute and store jackstraw plot
#' @param useAssay Assay containing scaled counts to use in JackStraw calculation.
#' @param dims Number of components to test in Jackstraw. If \code{NULL}, then all components are used. Default \code{NULL}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' sce <- seuratScaleData(sce, useAssay = "counts")
#' sce <- seuratPCA(sce, useAssay = "counts")
#' sce <- seuratComputeJackStraw(sce, useAssay = "counts")
#' }
#' @return Updated \code{SingleCellExperiment} object with jackstraw computations stored in it
#' @export
seuratComputeJackStraw <- function(inSCE, useAssay, dims = NULL) {
    seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
    if(is.null(seuratObject@reductions[["pca"]])) {
      stop("'seuratPCA' must be run before JackStraw can be computed.")
    }
    if(is.null(dims)) {
      dims <- ncol(seuratObject@reductions[["pca"]])
    }
    seuratObject <- Seurat::JackStraw(seuratObject, dims = as.double(dims))
    seuratObject <- Seurat::ScoreJackStraw(seuratObject, dims = 1:dims)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratJackStrawPlot
#' Computes the plot object for jackstraw plot from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute the jackstraw plot (pca should be computed)
#' @param dims Number of components to plot in Jackstraw. If \code{NULL}, then all components are plotted Default \code{NULL}.
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
seuratJackStrawPlot <- function(inSCE, dims = NULL) {
  seuratObject <- convertSCEToSeurat(inSCE)
  if(is.null(seuratObject@reductions[["pca"]])) {
    stop("'seuratPCA' must be run before JackStraw can be computed.")
  }
  if(is.null(dims)) {
    dims <- ncol(seuratObject@reductions[["pca"]])
  }
  return(Seurat::JackStrawPlot(seuratObject, dims = 1:dims))
}

#' seuratPlotHVG
#' Plot highly variable genes from input sce object (must have highly variable genes computations stored)
#' @param inSCE (sce) object that contains the highly variable genes computations
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- seuratNormalizeData(sce, useAssay = "counts")
#' sce <- seuratFindHVG(sce, useAssay = "counts")
#' seuratPlotHVG(sce)
#' }
#' @return plot object
#' @export
seuratPlotHVG <- function(inSCE) {
    seuratObject <- convertSCEToSeurat(inSCE)
    plot <- Seurat::VariableFeaturePlot(seuratObject)
    plot$labels$colour <- "Variable"
    if(requireNamespace("stringr", quietly = TRUE)){
      plot$data$colors <- stringr::str_to_title(plot$data$colors)
    }
    return(plot)
}

#' seuratReductionPlot
#' Plots the selected dimensionality reduction method
#' @param inSCE (sce) object which has the selected dimensionality reduction algorithm already computed and stored
#' @param useReduction Dimentionality reduction to plot. One of "pca", "ica", "tsne", or "umap". Default \code{"umap"}.
#' @param showLegend Select if legends should be shown on the output plot or not. Either "TRUE" or "FALSE". Default \code{FALSE}.
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
seuratReductionPlot <- function(inSCE, useReduction = c("pca", "ica", "tsne", "umap"), showLegend = FALSE) {
    seuratObject <- convertSCEToSeurat(inSCE)
    if(showLegend){
      plot <- Seurat::DimPlot(seuratObject, reduction = useReduction)
    }else{
      plot <- Seurat::DimPlot(seuratObject, reduction = useReduction) + Seurat::NoLegend()
    }
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}


#' seuratFindClusters
#' Computes the clusters from the input sce object and stores them back in sce object
#' @param inSCE (sce) object from which clusters should be computed and stored in
#' @param useAssay Assay containing scaled counts to use for clustering.
#' @param useReduction Reduction method to use for computing clusters. One of "pca" or "ica". Default \code{"pca"}.
#' @param dims numeric value of how many components to use for computing clusters. Default \code{10}.
#' @param algorithm selected algorithm to compute clusters. One of "louvain", "multilevel", or "SLM". Use \code{louvain} for "original Louvain algorithm" and \code{multilevel} for "Louvain algorithm with multilevel refinement". Default \code{louvain}.
#' @param groupSingletons boolean if singletons should be grouped together or not. Default \code{TRUE}.
#' @param resolution Set the resolution parameter to find larger (value above 1) or smaller (value below 1) number of communities. Default \code{0.8}.
#' @param externalReduction Pass DimReduc object if PCA/ICA computed through other libraries. Default \code{NULL}.
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
seuratFindClusters <- function(inSCE, useAssay, useReduction = c("pca", "ica"), dims = 10, algorithm = c("louvain", "multilevel", "SLM"), groupSingletons = TRUE, resolution = 0.8, externalReduction = NULL) {


    algorithm <- match.arg(algorithm)
    useReduction <- match.arg(useReduction)

    seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)

    if(!is.null(externalReduction)){
      seuratObject@reductions <- list(pca = externalReduction)
    }

    seuratObject <- Seurat::FindNeighbors(seuratObject, reduction = useReduction, dims = seq(dims))
    no_algorithm <- 1
    if (algorithm == "louvain") {
        no_algorithm = 1
    } else if (algorithm == "multilevel") {
        no_algorithm = 2
    } else if (algorithm == "SLM") {
        no_algorithm = 3
    }
    seuratObject <- Seurat::FindClusters(seuratObject, algorithm = no_algorithm, group.singletons = groupSingletons, resolution = resolution)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    colData(inSCE)[[paste0("Seurat","_",algorithm,"_","Resolution",resolution)]] <- seuratObject@meta.data$seurat_clusters
    return(inSCE)
}

#' seuratRunTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back into the sce object
#' @param inSCE (sce) object on which to compute the tSNE
#' @param useReduction selected reduction algorithm to use for computing tSNE. One of "pca" or "ica". Default \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Seurat tSNE Default \code{seuratTSNE}.
#' @param dims Number of reduction components to use for tSNE computation. Default \code{10}.
#' @param perplexity Adjust the preplexity tuneable parameter for the underlying tSNE call. Default \code{30}.
#' @return Updated sce object with tSNE computations stored
#' @export
#' @importFrom SingleCellExperiment reducedDim<-
seuratRunTSNE <- function(inSCE, useReduction = c("pca", "ica"), reducedDimName = "seuratTSNE", dims = 10, perplexity = 30) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE)
  seuratObject <- Seurat::RunTSNE(seuratObject, reduction = useReduction, dims = 1:dims, perplexity = perplexity)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  temp <- seuratObject@reductions$tsne@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#RunUMAP(seurat, reduction = "pca", dims = 1:10, min.dist = 0.4, n.neighbors = 40, spread = 20)
#' seuratRunUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back into the sce object
#' @param inSCE (sce) object on which to compute the UMAP
#' @param useReduction Reduction to use for computing UMAP. One of "pca" or "ica". Default is \code{"pca"}.
#' @param reducedDimName Name of new reducedDims object containing Seurat UMAP Default \code{seuratUMAP}.
#' @param dims Numerical value of how many reduction components to use for UMAP computation. Default \code{10}.
#' @param minDist Sets the \code{"min.dist"} parameter to the underlying UMAP call. See \link[Seurat]{RunUMAP} for more information. Default \code{0.3}.
#' @param nNeighbors Sets the \code{"n.neighbors"} parameter to the underlying UMAP call. See \link[Seurat]{RunUMAP} for more information. Default \code{30L}.
#' @param spread Sets the \code{"spread"} parameter to the underlying UMAP call. See \link[Seurat]{RunUMAP} for more information. Default \code{1}.
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
seuratRunUMAP <- function(inSCE, useReduction = c("pca", "ica"), reducedDimName = "seuratUMAP", dims = 10, minDist = 0.3, nNeighbors = 30L, spread = 1) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE)
  seuratObject <- Seurat::RunUMAP(seuratObject,
                                  reduction = useReduction,
                                  dims = 1:dims,
                                  min.dist = minDist,
                                  n.neighbors = nNeighbors,
                                  spread = spread)
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)

  temp <- seuratObject@reductions$umap@cell.embeddings
  rownames(temp) <- colnames(inSCE)
  reducedDim(inSCE, reducedDimName) <- temp

  return(inSCE)
}

#' .seuratGetVariableFeatures
#' Retrieves the requested number of variable feature names
#' @param inSCE (sce) object from which to extract the variable feature names
#' @param numberOfFeatures numerical value indicating how many feature names should be retrieved (default is 100)
#' @return list() of variable feature names
.seuratGetVariableFeatures <- function(inSCE, numberOfFeatures) {
    seuratObject <- convertSCEToSeurat(inSCE)
    if (length(seuratObject@assays$RNA@var.features) > 0) {
        return(seuratObject@assays$RNA@var.features[1:numberOfFeatures])
    }
}

#' seuratElbowPlot
#' Computes the plot object for elbow plot from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute the elbow plot (pca should be computed)
#' @param significantPC Number of significant principal components to plot. This is used to alter the color of the points for the corresponding PCs. If \code{NULL}, all points will be the same color. Default \code{NULL}.
#' @param reduction Reduction to use for elbow plot generation. Either \code{"pca"} or \code{"ica"}. Default \code{"pca"}.
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
seuratElbowPlot <- function(inSCE, significantPC = NULL, reduction = "pca") {
    seuratObject <- convertSCEToSeurat(inSCE)
    plot <- Seurat::ElbowPlot(seuratObject, reduction = reduction)
    if(!is.null(significantPC)){
      plot$data$Significant <- c(rep("Yes", significantPC), rep("No", length(rownames(plot$data)) - significantPC))
      plot <- ggplot2::ggplot(data = plot$data, ggplot2::aes(x = plot$data$dims, y = plot$data$stdev, color = plot$data$Significant)) + ggplot2::geom_point()
    }
    plot$labels$x <- "PC"
    plot$labels$y <- "Standard Deviation"
    plot$labels$colour <- "Significant"
    return(plot)
}

#' seuratComputeHeatmap
#' Computes the heatmap plot object from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute heatmap (pca should be computed)
#' @param useAssay Assay containing scaled counts to use in heatmap.
#' @param useReduction Reduction method to use for computing clusters. One of "pca" or "ica". Default \code{"pca"}.
#' @param dims Number of components to generate heatmap plot objects. If \code{NULL}, a heatmap will be generated for all components. Default \code{NULL}.
#' @param nfeatures Numer of features to include in the heatmap. Default \code{30}.
#' @param fast See \link[Seurat]{DimHeatmap} for more information. Default \code{TRUE}.
#' @param combine See \link[Seurat]{DimHeatmap} for more information. Default \code{TRUE}.
#' @param raster See \link[Seurat]{DimHeatmap} for more information. Default \code{TRUE}.
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
seuratComputeHeatmap <- function(inSCE, useAssay, useReduction = c("pca", "ica"), dims = NULL, nfeatures = 30, fast = TRUE, combine = TRUE, raster = TRUE) {
  useReduction <- match.arg(useReduction)
  seuratObject <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  if(is.null(dims)) {
    dims <- ncol(seuratObject@reductions[[useReduction]])
  }
  return(Seurat::DimHeatmap(seuratObject, dims = 1:dims, nfeatures = nfeatures, reduction = useReduction, fast = fast, combine = combine, raster = raster))
}

#' seuratHeatmapPlot
#' Modifies the heatmap plot object so it contains specified number of heatmaps in a single plot
#' @param plotObject plot object computed from seuratComputeHeatmap() function
#' @param dims numerical value of how many heatmaps to draw (default is 0)
#' @param ncol numerical value indicating that in how many columns should the heatmaps be distrbuted (default is 2)
#' @param labels list() of labels to draw on heatmaps
#' @return modified plot object
#' @export
seuratHeatmapPlot <- function(plotObject, dims, ncol, labels) {
    componentsToPlot <- as.integer(gsub("[^0-9.]", "", labels))
    return(cowplot::plot_grid(plotlist = plotObject[c(componentsToPlot)], ncol = ncol, labels = labels))
}

#' .updateAssaySCE
#' Update/Modify/Add an assay in the provided SingleCellExperiment object from a Seurat object
#' @param inSCE Input SingleCellExperiment object
#' @param seuratObject Input Seurat object
#' @param assaySlotSCE Selected assay to update in the input SingleCellExperiment object
#' @param seuratDataSlot Selected data slot from the Seurat object. Default \code{"counts"}.
#' @param seuratAssaySlot Selected assay from Seurat object. Default \code{"RNA"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  data from Seurat object appended to the \link[SummarizedExperiment]{assay} slot.
#' @importFrom SummarizedExperiment assay<-
.updateAssaySCE <- function(inSCE, seuratObject, assaySlotSCE, seuratDataSlot = "counts", seuratAssaySlot = "RNA") {
  assay(inSCE, assaySlotSCE) <- NULL
  temp.matrix <- methods::slot(Seurat::GetAssay(seuratObject, seuratAssaySlot), seuratDataSlot)
  rownames(temp.matrix) <- rownames(inSCE)
  colnames(temp.matrix) <- colnames(inSCE)
  assay(inSCE, assaySlotSCE) <- temp.matrix
  return(inSCE)
}

#' convertSeuratToSCE
#' Converts the input seurat object to a sce object
#' @param seuratObject Input Seurat object
#' @param normAssayName Name of assay to store the normalized data. Default \code{"seuratNormData"}.
#' @param scaledAssayName Name of assay to store the scaled data. Default \code{"seuratScaledData"}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' seurat <- convertSCEToSeurat(sce)
#' sce <- convertSeuratToSCE(seurat)
#' @return \code{SingleCellExperiment} output object
#' @export
convertSeuratToSCE <- function(seuratObject, normAssayName = "seuratNormData", scaledAssayName = "seuratScaledData") {
  inSCE <- Seurat::as.SingleCellExperiment(seuratObject)
  assay(inSCE, normAssayName) <- methods::slot(seuratObject@assays$RNA, "data")
  if (length(methods::slot(seuratObject, "assays")[["RNA"]]@scale.data) > 0) {
    assay(inSCE, scaledAssayName) <- methods::slot(seuratObject@assays$RNA, "scale.data")
  }
  inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
  return(inSCE)
}

#' convertSCEToSeurat
#' Converts sce object to seurat while retaining all assays and metadata
#' @param inSCE A \code{SingleCellExperiment} object to convert to a Seurat object.
#' @param countsAssay Which assay to use from sce object for raw counts. Default \code{NULL}.
#' @param normAssay Which assay to use from sce object for normalized data. Default \code{NULL}.
#' @param scaledAssay Which assay to use from sce object for scaled data. Default \code{NULL}
#' @examples
#' data(scExample, package = "singleCellTK")
#' seurat <- convertSCEToSeurat(sce)
#' @return Updated seurat object that contains all data from the input sce object
#' @export
#' @importFrom SummarizedExperiment assay assays
convertSCEToSeurat <- function(inSCE, countsAssay = NULL, normAssay = NULL, scaledAssay = NULL) {

  if(!is.null(countsAssay) && !(countsAssay %in% names(assays(inSCE)))) {
    stop(paste0("'", countsAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
  }
  if(!is.null(normAssay) && !(normAssay %in% names(assays(inSCE)))) {
    stop(paste0("'", normAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
  }
  if(!is.null(scaledAssay) && !(scaledAssay %in% names(assays(inSCE)))) {
    stop(paste0("'", scaledAssay, "' not found in the list of assays: ",
                paste(names(assays(inSCE)), collapse=",")))
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
    seuratObject@assays$RNA@data <- .convertToMatrix(assay(inSCE, normAssay))
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
  return(seuratObject)
}

#' seuratSCTransform
#' Runs the \link[Seurat]{SCTransform} function to transform/normalize the input data
#' @param inSCE Input SingleCellExperiment object
#' @param normAssayName Name for the output data assay. Default \code{"SCTCounts"}.
#' @param useAssay Name for the input data assay. Default \code{"counts"}.
#' @return Updated SingleCellExperiment object containing the transformed data
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- seuratSCTransform(sce_chcl, "SCTCounts", "counts")
seuratSCTransform <- function(inSCE, normAssayName = "SCTCounts", useAssay = "counts") {
    seuratObject <- base::suppressWarnings(Seurat::SCTransform(
                            object = convertSCEToSeurat(inSCE, useAssay),
                            assay = "RNA",
                            new.assay.name = "SCTransform",
                            do.correct.umi = FALSE,
                            verbose = TRUE))
    inSCE <- .updateAssaySCE(inSCE = inSCE, seuratObject = seuratObject, assaySlotSCE = normAssayName, seuratDataSlot = "data", seuratAssaySlot = "SCTransform")
    return(inSCE)
}


#' .seuratInvalidate
#' Removes seurat data from the input SingleCellExperiment object specified by the task in the Seurat workflow.
#' @param inSCE Input \code{SingleCellExperiment} object to remove Seurat data from.
#' @param scaleData Remove scaled data from seurat. Default \code{TRUE}.
#' @param varFeatures Remove variable features from seurat. Default \code{TRUE}.
#' @param PCA Remove PCA from seurat. Default \code{TRUE}.
#' @param ICA Remove ICA from seurat. Default \code{TRUE}.
#' @param tSNE Remove tSNE from seurat. Default \code{TRUE}.
#' @param UMAP Remove UMAP from seurat. Default \code{TRUE}.
#' @param clusters Remove clusters from seurat. Default \code{TRUE}.
#' @return Updated SingleCellExperiment object containing the Seurat object in the metadata slot with the data removed
#' @importFrom SummarizedExperiment assay<-
.seuratInvalidate <- function(inSCE, scaleData = TRUE, varFeatures = TRUE, PCA = TRUE, ICA = TRUE, tSNE = TRUE, UMAP = TRUE, clusters = TRUE){
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
#' @param inSCE Input \code{SingleCellExperiment} object that contains the assay to batch-correct.
#' @param useAssay Assay to batch-correct.
#' @param batch Batch variable from \code{colData} slot of \code{SingleCellExperiment} object.
#' @param newAssayName Assay name for the batch-corrected ouput assay.
#' @param kAnchor Number of neighbours to use for finding the anchors in the \link[Seurat]{FindIntegrationAnchors} function.
#' @param kFilter Number of neighbours to use for filtering the anchors in the \link[Seurat]{FindIntegrationAnchors} function.
#' @param kWeight Number of neighbours to use when weigthing the anchors in the \link[Seurat]{IntegrateData} function.
#' @param ndims Number of dimensions to use. Default \code{10}.
#'
#' @return A \code{SingleCellExperiment} object that contains the batch-corrected assay inside the \code{altExp} slot of the object
#' @export
seuratIntegration <- function(inSCE, useAssay = "counts", batch, newAssayName = "SeuratIntegratedAssay", kAnchor, kFilter, kWeight, ndims = 10){
   if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
    stop(paste(useAssay, "not found in the input object assays"))
  }
  if(is.null(batch)){
    stop("batch variable must be provided for batch-correction")
  }
  if(kAnchor == 0 || kFilter == 0 || kWeight == 0){
    stop("kAnchor, kFilter or kWeight cannot be zero. Please input correct parameters.")
  }

  #create seurat object
  seuratObject <- convertSCEToSeurat(inSCE, useAssay)
  rownames(seuratObject@meta.data) <- gsub("_", "-", rownames(seuratObject@meta.data))

  #split seurat object by batch variable
  seurat.list <- Seurat::SplitObject(seuratObject, split.by = batch)
  seurat.list <- seurat.list[c(unique(seuratObject@meta.data[[batch]]))]

  #find anchors
  seurat.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat.list, dims = 1:ndims, k.anchor = kAnchor, k.filter = kFilter)
  seurat.integrated <- Seurat::IntegrateData(anchorset = seurat.anchors, dims = 1:ndims, k.weight = kWeight)

  #store results back in altExp slot of sce object
  altExp(inSCE, newAssayName) <- SingleCellExperiment(list(counts = Seurat::GetAssayData(seurat.integrated@assays$integrated, "data")))
  SummarizedExperiment::assayNames(altExp(inSCE,newAssayName)) <- newAssayName #remove this if counts in above line set to altExp

  #store back colData from sce into the altExp slot
  colData(altExp(inSCE, newAssayName))<- colData(inSCE)

  #counts <- assay(altExp(inSCE, newAssayName), "altExp")
  #remove NA values from counts and replace with zero so can be used properly by dgCMatrix
  counts <- assay(altExp(inSCE, newAssayName), newAssayName)
  counts[is.na(counts)] <- 0

  #store back counts
  #assay(altExp(inSCE, newAssayName), "altExp") <- counts
  assay(altExp(inSCE, newAssayName), newAssayName) <- counts

  return(inSCE)
}

