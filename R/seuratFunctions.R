# Helper/Wrapper Functions ---

#' .getPCAComponentNames
#' Creates a list of PC components to populate the picker for PC heatmap generation
#' @param maxComponents; number of components to return for the picker
#' @return componentNames; list of component names (appended with "PC")
.getPCAComponentNames <- function(maxComponents) {
    componentNames <- list()
    for (i in 1:maxComponents) {
        componentNames[i] <- paste0("PC", i)
    }
    return(componentNames)
}

#' .rdsToSce
#' Reads rds file (from a local path) and loads into sce object 
#' *Only to be used for first time initialization of the rds file into R environment*
#' @param filePath; path of the rds file to load
#' @return sce object
.rdsToSce <- function(filePath) {
    sce <- readRDS(filePath)
    return(sce)
}

#' .sceToSeurat
#' Converts a sce object to seurat object (using rds filepath)
#' *Only to be used for first time initialization of seurat object*
#' @param filePath; path of the rds file to convert to seurat object
#' @return seurat object
.sceToSeurat <- function(filePath) {
    seuratObject <- CreateSeuratObject(counts = counts(.rdsToSce(filePath)))
    return(seuratObject)
}

#' .addSeuratToMetaDataSCE
#' Adds the input seurat object to the metadata slot of the input sce object (after removing the data matrices)
#' @param sce; sce object to which seurat object should be added in the metadata slot (copy to)
#' @param seuratObject; seurat object which should be added to the metadata slot of sce object (copy from)
#' @return sce; updated sce object which now contains the seurat object in its metadata slot (excluding data matrices)
.addSeuratToMetaDataSCE <- function(sce, seuratObject) {
    seuratObject@assays$RNA@counts <- new("dgCMatrix")
    seuratObject@assays$RNA@data <- new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    sce@metadata[["seurat"]] <- seuratObject
    return(sce)
}

#' .rowNamesSeurat
#' Retrieves a list of genenames/rownames/featurenames from seurat object
#' @param seuratObject; seurat object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSeurat <- function(seuratObject) {
    return(rownames(seuratObject))
}

#' .rowNamesSCE
#' Retrieves a list of genenames/rownames/featurenames from sce object
#' @param sce; sce object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSCE <- function(sce) {
    return(rownames(sce))
}

#' .computeSignificantPC
#' Computes the significant principal components from an input sce object (must containt pca slot) using stdev
#' @param sceObject; sce object with pca computed
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return max_components; a numerical value indicating how many number of components are considered significant
.computeSignificantPC <- function(sceObject, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    max_components <- 0
    for (i in 1:(length(seuratObject[["pca"]]@stdev) - 1)) {
        if (abs(seuratObject[["pca"]]@stdev[i + 1] - seuratObject[["pca"]]@stdev[i]) > 0.1) {
            max_components <- i
        }
    }
    return(max_components)
}

#' seuratNormalizeData
#' Wrapper for NormalizeData() function from seurat library
#' Normalizes the sce object according to the input parameters 
#' @param sceObject; sce object to normalize
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param normalizationMethod; selected normalization method (default is "LogNormalize")
#' @param scaleFactor; numeric value that represents the scaling factor (default is 10000)
#' @return sceObject; normalized sce object
#' @export
seuratNormalizeData <- function(sceObject, geneNamesSeurat, normalizationMethod, scaleFactor) {
    seuratObject <- NormalizeData(convertSCEToSeurat(sceObject, geneNamesSeurat), normalization.method = normalizationMethod, scale.factor = scaleFactor)
    sceObject <- convertSeuratToSCE(sceObject, geneNamesSeurat, seuratObject, "seuratNormalizedData", "data")
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratScaleData
#' Scales the input sce object according to the input parameters
#' @param sceObject; sce object to scale
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param model.use; selected model to use for scaling data (default is "linear")
#' @param do.scale; boolean if data should be scaled or not (TRUE or FALSE, default is TRUE)
#' @param do.center; boolean if data should be centered or not (TRUE or FALSE, default is TRUE)
#' @param scale.max; maximum numeric value to return for scaled data (default is 10)
#' @return sceObject; scaled sce object
#' @export
seuratScaleData <- function(sceObject, geneNamesSeurat, model.use, do.scale, do.center, scale.max) {
    seuratObject <- ScaleData(convertSCEToSeurat(sceObject, geneNamesSeurat), model.use = model.use, do.scale = do.scale, do.center = do.center, scale.max = as.double(scale.max))
    sceObject <- convertSeuratToSCE(sceObject, geneNamesSeurat, seuratObject, "seuratScaledData", "scale.data")
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratPCA
#' Computes PCA on the input sce object and stores the calculated principal components within the sce object
#' @param sceObject; sce object on which to compute PCA
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param npcs; numeric value of how many components to compute (default is 20)
#' @return sceObject; updated sce object which now contains the computed principal components
#' @export
seuratPCA <- function(sceObject, geneNamesSeurat, npcs) {
    seuratObject <- RunPCA(convertSCEToSeurat(sceObject, geneNamesSeurat), npcs = as.double(npcs))
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratICA
#' Computes ICA on the input sce object and stores the calculated independent components within the sce object
#' @param sceObject; sce object on which to compute ICA
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param nics; numeric value of how many components to compute (default is 20)
#' @return sceObject; updated sce object which now contains the computed independent components
#' @export
seuratICA <- function(sceObject, geneNamesSeurat, nics) {
    seuratObject <- RunICA(convertSCEToSeurat(sceObject, geneNamesSeurat), nics = as.double(nics))
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratComputeJackStraw
#' Compute jackstraw plot and store the computations in the input sce object
#' @param sceObject; sce object on which to compute and store jackstraw plot
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; numeric value of how many components to use for jackstraw plot (default = number of computed principal components)
#' @return sceObject; updated sce object with jackstraw computations stored in it
#' @export
seuratComputeJackStraw <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- JackStraw(seuratObject, dims = as.double(dims))
    seuratObject <- ScoreJackStraw(seuratObject, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratFindHVG
#' Find highly variable genes and store in the input sce object
#' @param sceObject; sce object to compute highly variable genes from and to store back to it
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param hvgMethod; selected method to use for computation of highly variable genes (default is "vst")
#' @param hvgNumber; numeric value of how many genes to select as highly variable (default is 2000)
#' @return sceObject; updated sce object with highly variable genes computation stored
#' @export
seuratFindHVG <- function(sceObject, geneNamesSeurat, hvgMethod, hvgNumber) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratPlotHVG
#' Plot highly variable genes from input sce object (must have highly variable genes computations stored)
#' @param sceObject; sce object that contains the highly variable genes computations
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return plot object 
#' @export
seuratPlotHVG <- function(sceObject, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(VariableFeaturePlot(seuratObject))
}

#' seuratReductionPlot
#' Plots the selected dimensionality reduction method
#' @param sceObject; sce object which has the selected dimensionality reduction algorithm already computed and stored
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; one of selected algorithm from pca, ica, tsne and umap
#' @return plot object
#' @export
seuratReductionPlot <- function(sceObject, geneNamesSeurat, reduction) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    plot <- DimPlot(seuratObject, reduction = reduction)
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}


#' convertSeuratToSCE
#' Modifies the input sce object to include the updated assays from seurat object
#' @param sce; outdated sce object in which we wish to include the newly calculated assays from seurat object
#' @param geneNames; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param seuratObject; updated seurat object that contains the newly calculated assays that we wish to store in the sce object
#' @param assaySlotSCE; the relevant assay slot in sce object (copy to)
#' @param assaySlotSeurat; the relevant assay slow in seurat object (copy from)
#' @return sce; sce object that contains the newly added/modified assays from the seurat object
#' @export
convertSeuratToSCE <- function(sce, geneNames, seuratObject, assaySlotSCE, assaySlotSeurat) {
    assay(sce, assaySlotSCE) <- NULL
    assay(sce, assaySlotSCE) <- slot(seuratObject@assays$RNA, assaySlotSeurat)
    rownames(sce) <- geneNames
    return(sce)
}

#' convertSCEToSeurat
#' Converts sce object to seurat while retaining all assays and metadata
#' @param sce; sce object to convert to seurat
#' @param geneNames; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return seuratObject; updated seurat object that contains all data from the input sce object
#' @export
convertSCEToSeurat <- function(sce, geneNames) {
    seuratObject <- CreateSeuratObject(counts = counts(sce))
    if ("seuratNormalizedData" %in% names(assays(sce))) {
        seuratObject@assays$RNA@data <- assay(sce, "seuratNormalizedData")
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if ("seuratScaledData" %in% names(assays(sce))) {
        seuratObject@assays$RNA@scale.data <- assay(sce, "seuratScaledData")
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if (!is.null(sce@metadata[["seurat"]]) && length(sce@metadata[["seurat"]]@assays$RNA@var.features) > 0) {
        seuratObject@assays$RNA@var.features <- sce@metadata[["seurat"]]@assays$RNA@var.features
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$pca)) {
        seuratObject@reductions$pca <- sce@metadata[["seurat"]]@reductions$pca
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@assays$RNA@meta.features)) {
        seuratObject@assays$RNA@meta.features <- sce@metadata[["seurat"]]@assays$RNA@meta.features
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$ica)) {
        seuratObject@reductions$ica <- sce@metadata[["seurat"]]@reductions$ica
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$tsne)) {
        seuratObject@reductions$tsne <- sce@metadata[["seurat"]]@reductions$tsne
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$umap)) {
        seuratObject@reductions$umap <- sce@metadata[["seurat"]]@reductions$umap
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@meta.data)) {
        seuratObject@meta.data <- sce@metadata[["seurat"]]@meta.data
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@commands)) {
        seuratObject@commands <- sce@metadata[["seurat"]]@commands
    }
    return(seuratObject)
}

#' seuratFindClusters
#' Computes the clusters from the input sce object and stores them back in sce object
#' @param sceObject; sce object from which clusters should be computed and stored in
#' @param reduction; selected reduction method to use for computing clusters ("pca" or "ica", default is "pca")
#' @param dims; numeric value of how many components to use for computing clusters (default is 10)
#' @param algorithm; selected algorithm to compute clusters (default is "original Louvain algorithm")
#' @param group.singletons; boolean if singletons should be grouped together or not (TRUE or FALSE, default is TRUE)
#' @return sceObject; updated sce object which now contains the computed clusters
#' @export
seuratFindClusters <- function(sceObject, geneNamesSeurat, reduction, dims, algorithm, group.singletons) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:dims)
    no_algorithm <- 1
    if (algorithm == "original Louvain algorithm") {
        no_algorithm = 1
    } else if (algorithm == "Louvain algorithm with multilevel refinement") {
        no_algorithm = 2
    } else if (algorithm == "SLM algorithm") {
        no_algorithm = 3
    }
    seuratObject <- FindClusters(seuratObject, algorithm = no_algorithm, group.singletons = group.singletons)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratRunTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back into the sce object
#' @param sceObject; sce object on which to compute the tSNE
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; selected reduction algorithm to use for computing tSNE ("pca" or "ica", default is "pca")
#' @param dims; numerical value of how many reduction components to use for tSNE computation (default is 10)
#' @return sceObject; updated sce object with tSNE computations stored
#' @export
seuratRunTSNE <- function(sceObject, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- RunTSNE(seuratObject, reduction = reduction, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratRunUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back into the sce object
#' @param sceObject; sce object on which to compute the UMAP
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; selected reduction algorithm to use for computing UMAP ("pca" or "ica", default is "pca")
#' @param dims; numerical value of how many reduction components to use for UMAP computation (default is 10)
#' @return sceObject; updated sce object with UMAP computations stored
#' @export
seuratRunUMAP <- function(sceObject, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- RunUMAP(seuratObject, reduction = reduction, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' .seuratGetVariableFeatures
#' Retrieves the requested number of variable feature names
#' @param sceObject; sce object from which to extract the variable feature names
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param numberOfFeatures; numerical value indicating how many feature names should be retrieved (default is 100)
#' @return list() of variable feature names
.seuratGetVariableFeatures <- function(sceObject, geneNamesSeurat, numberOfFeatures) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    if (length(seuratObject@assays$RNA@var.features) > 0) {
        return(print(seuratObject@assays$RNA@var.features[1:numberOfFeatures]))
    }
}

#' seuratElbowPlot
#' Computes the plot object for elbow plot from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute the elbow plot (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param significantPC; a numerical value indicating the number of significant principal components (used to alter the color of the significant components)
#' @return plot object
#' @export
seuratElbowPlot <- function(sceObject, geneNamesSeurat, significantPC) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    plot <- ElbowPlot(seuratObject)
    plot <- ggplot_build(plot)
    for (i in 1:significantPC) {
        plot$data[[1]]$shape[i] <- 16
        plot$data[[1]]$colour[i] <- "red"
        plot$data[[1]]$size[i] <- 3.5
    }
    plot <- ggplot_gtable(plot)
    plot <- as.ggplot(plot)
    return(plot)
}

#' seuratJackStrawPlot
#' Computes the plot object for jackstraw plot from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute the jackstraw plot (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; a numerical value indicating how many components to use in the computation of jackstraw plot from pca (default is number of pca components computed)
#' @return plot object
#' @export
seuratJackStrawPlot <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(JackStrawPlot(seuratObject, dims = 1:dims))
}

#' seuratComputeHeatmap
#' Computes the heatmap plot object from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute heatmap (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; numerical value indicating how many components to use for the computation of heatmap plot object (default is number of pca components computed) 
#' @return plot object
#' @export
seuratComputeHeatmap <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(DimHeatmap(seuratObject, dims = 1:dims, fast = FALSE, combine = FALSE))
}

#' seuratHeatmapPlot
#' Modifies the heatmap plot object so it contains specified number of heatmaps in a single plot 
#' @param plotObject; plot object computed from seuratComputeHeatmap() function
#' @param dims; numerical value of how many heatmaps to draw (default is 0)
#' @param ncol; numerical value indicating that in how many columns should the heatmaps be distrbuted (default is 2)
#' @param labels; list() of labels to draw on heatmaps
#' @return modified plot object
#' @export
seuratHeatmapPlot <- function(plotObject, dims, ncol, labels) {
    componentsToPlot <- as.integer(gsub("[^0-9.]", "", labels))
    return(plot_grid(plotlist = plotObject[c(componentsToPlot)], ncol = ncol, labels = labels))
}

# ----