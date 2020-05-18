# Helper/Wrapper Functions ---

#' .getPCAComponentNames
#' Creates a list of PC components to populate the picker for PC heatmap generation
#' @param maxComponents number of components to return for the picker
#' @return componentNames list of component names (appended with "PC")
.getPCAComponentNames <- function(maxComponents) {
    componentNames <- list()
    for (i in seq(maxComponents)) {
        componentNames[i] <- paste0("PC", i)
    }
    return(componentNames)
}

#' .rdsToSce
#' Reads rds file (from a local path) and loads into sce object 
#' *Only to be used for first time initialization of the rds file into R environment*
#' @param filePath path of the rds file to load
#' @return sce object
.rdsToSce <- function(filePath) {
    inSCE <- readRDS(filePath)
    return(inSCE)
}

#' .sceToSeurat
#' Converts a sce object to seurat object (using rds filepath)
#' *Only to be used for first time initialization of seurat object*
#' @param filePath path of the rds file to convert to seurat object
#' @return seurat object
.sceToSeurat <- function(filePath) {
    seuratObject <- Seurat::CreateSeuratObject(counts = counts(.rdsToSce(filePath)))
    return(seuratObject)
}

#' .addSeuratToMetaDataSCE
#' Adds the input seurat object to the metadata slot of the input sce object (after removing the data matrices)
#' @param inSCE (sce) object to which seurat object should be added in the metadata slot (copy to)
#' @param seuratObject seurat object which should be added to the metadata slot of sce object (copy from)
#' @return updated sce object which now contains the seurat object in its metadata slot (excluding data matrices)
.addSeuratToMetaDataSCE <- function(inSCE, seuratObject) {
    seuratObject@assays$RNA@counts <- methods::new("dgCMatrix")
    seuratObject@assays$RNA@data <- methods::new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    inSCE@metadata[["seurat"]] <- seuratObject
    return(inSCE)
}

#' .rowNamesSeurat
#' Retrieves a list of genenames/rownames/featurenames from seurat object
#' @param seuratObject seurat object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSeurat <- function(seuratObject) {
    return(rownames(seuratObject))
}

#' .rowNamesSCE
#' Retrieves a list of genenames/rownames/featurenames from sce object
#' @param inSCE sce object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSCE <- function(inSCE) {
    return(rownames(inSCE))
}

#' .computeSignificantPC
#' Computes the significant principal components from an input sce object (must containt pca slot) using stdev
#' @param inSCE (sce) object with pca computed
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return max_components a numerical value indicating how many number of components are considered significant
.computeSignificantPC <- function(inSCE, useAssay, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
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
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param normalizationMethod selected normalization method (default is "LogNormalize")
#' @param scaleFactor numeric value that represents the scaling factor (default is 10000)
#' @return sceObject normalized sce object
#' @export
seuratNormalizeData <- function(inSCE, newAssayName = "seuratNormalizedData", useAssay, geneNamesSeurat, normalizationMethod = "LogNormalize", scaleFactor = 10000) {
    seuratObject <- Seurat::NormalizeData(convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat), normalization.method = normalizationMethod, scale.factor = scaleFactor)
    inSCE <- .updateAssaySCE(inSCE, geneNamesSeurat, seuratObject, newAssayName, "data")
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    inSCE@metadata$selected_assay <- useAssay
    return(inSCE)
}

#' seuratScaleData
#' Scales the input sce object according to the input parameters
#' @param inSCE (sce) object to scale
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param model.use selected model to use for scaling data (default is "linear")
#' @param do.scale boolean if data should be scaled or not (TRUE or FALSE, default is TRUE)
#' @param do.center boolean if data should be centered or not (TRUE or FALSE, default is TRUE)
#' @param scale.max maximum numeric value to return for scaled data (default is 10)
#' @return sceObject scaled sce object
#' @export
seuratScaleData <- function(inSCE, useAssay, geneNamesSeurat, model.use = "linear", do.scale = TRUE, do.center = TRUE, scale.max = 10) {
    seuratObject <- Seurat::ScaleData(convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat), model.use = model.use, do.scale = do.scale, do.center = do.center, scale.max = as.double(scale.max))
    inSCE <- .updateAssaySCE(inSCE, geneNamesSeurat, seuratObject, "seuratScaledData", "scale.data")
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratPCA
#' Computes PCA on the input sce object and stores the calculated principal components within the sce object
#' @param inSCE (sce) object on which to compute PCA
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param nPCs numeric value of how many components to compute (default is 20)
#' @return sceObject updated sce object which now contains the computed principal components
#' @export
seuratPCA <- function(inSCE, useAssay, geneNamesSeurat, nPCs) {
    seuratObject <- Seurat::RunPCA(convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat), npcs = as.double(nPCs))
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratICA
#' Computes ICA on the input sce object and stores the calculated independent components within the sce object
#' @param inSCE (sce) object on which to compute ICA
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param nics numeric value of how many components to compute (default is 20)
#' @return sceObject updated sce object which now contains the computed independent components
#' @export
seuratICA <- function(inSCE, useAssay, geneNamesSeurat, nics) {
    seuratObject <- Seurat::RunICA(convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat), nics = as.double(nics))
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratComputeJackStraw
#' Compute jackstraw plot and store the computations in the input sce object
#' @param inSCE (sce) object on which to compute and store jackstraw plot
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims numeric value of how many components to use for jackstraw plot (default = number of computed principal components)
#' @return sceObject updated sce object with jackstraw computations stored in it
#' @export
seuratComputeJackStraw <- function(inSCE, useAssay, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    seuratObject <- Seurat::JackStraw(seuratObject, dims = as.double(dims))
    seuratObject <- Seurat::ScoreJackStraw(seuratObject, dims = 1:dims)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratFindHVG
#' Find highly variable genes and store in the input sce object
#' @param inSCE (sce) object to compute highly variable genes from and to store back to it
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param hvgMethod selected method to use for computation of highly variable genes (default is "vst")
#' @param hvgNumber numeric value of how many genes to select as highly variable (default is 2000)
#' @return sceObject updated sce object with highly variable genes computation stored
#' @export
seuratFindHVG <- function(inSCE, useAssay, geneNamesSeurat, hvgMethod, hvgNumber) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    seuratObject <- Seurat::FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    if (hvgMethod == "vst") {
        rowData(inSCE)$seurat_variableFeatures_vst_varianceStandardized <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$vst.variance.standardized
        rowData(inSCE)$seurat_variableFeatures_vst_mean <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$vst.mean
    } else if (hvgMethod == "dispersion") {
        rowData(inSCE)$seurat_variableFeatures_dispersion_dispersion <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.dispersion
        rowData(inSCE)$seurat_variableFeatures_dispersion_dispersionScaled <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.dispersion.scaled
        rowData(inSCE)$seurat_variableFeatures_dispersion_mean <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.mean
    }
    else if (hvgMethod == "mean.var.plot") {
        rowData(inSCE)$seurat_variableFeatures_mvp_dispersion <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.dispersion
        rowData(inSCE)$seurat_variableFeatures_mvp_dispersionScaled <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.dispersion.scaled
        rowData(inSCE)$seurat_variableFeatures_mvp_mean <- methods::slot(inSCE@metadata[["seurat"]], "assays")[["RNA"]]@meta.features$mvp.mean
    }
    return(inSCE)
}

#' seuratPlotHVG
#' Plot highly variable genes from input sce object (must have highly variable genes computations stored)
#' @param inSCE (sce) object that contains the highly variable genes computations
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return plot object 
#' @export
seuratPlotHVG <- function(inSCE, useAssay, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    return(Seurat::VariableFeaturePlot(seuratObject))
}

#' seuratReductionPlot
#' Plots the selected dimensionality reduction method
#' @param inSCE (sce) object which has the selected dimensionality reduction algorithm already computed and stored
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction one of selected algorithm from pca, ica, tsne and umap
#' @return plot object
#' @export
seuratReductionPlot <- function(inSCE, useAssay, geneNamesSeurat, reduction) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    plot <- Seurat::DimPlot(seuratObject, reduction = reduction)
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}

#' .updateAssaySCE
#' Update/Modify/Add an assay in the provided sce object
#' @param inSCE (sce) object to which we have to add assay to (copy to)
#' @param geneNames for consistent formatting
#' @param seuratObject from which we have to copy the assay (copy from)
#' @param assaySlotSCE the counts slot from assay in sce object
#' @param assaySlotSeurat the counts slot from given assay of seurat object
#' @param slotSeurat which assay to use in seurat object (by default it is "RNA")
.updateAssaySCE <- function(inSCE, geneNames, seuratObject, assaySlotSCE, assaySlotSeurat, slotSeurat = "RNA") {
    assay(inSCE, assaySlotSCE) <- NULL
    assay(inSCE, assaySlotSCE) <- methods::slot(Seurat::GetAssay(seuratObject, slotSeurat), assaySlotSeurat)
    rownames(inSCE) <- geneNames
    return(inSCE)
}

#' convertSeuratToSCE
#' Converts the input seurat object to a sce object
#' @param seuratObject input object
#' @return inSCE output object
#' @export
convertSeuratToSCE <- function(seuratObject) {
    inSCE <- Seurat::as.SingleCellExperiment(seuratObject)
    assay(inSCE, "seuratNormalizedData") <- methods::slot(seuratObject@assays$RNA, "data")
    if (length(methods::slot(seuratObject, "assays")[["RNA"]]@scale.data) > 0) {
        assay(inSCE, "seuratScaledData") <- methods::slot(seuratObject@assays$RNA, "scale.data")
    }
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' convertSCEToSeurat
#' Converts sce object to seurat while retaining all assays and metadata
#' @param inSCE (sce) object to convert to seurat
#' @param useAssay which assay to use from sce object
#' @param geneNames a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return seuratObject updated seurat object that contains all data from the input sce object
#' @export
convertSCEToSeurat <- function(inSCE, useAssay, geneNames) {
    seuratObject <- Seurat::CreateSeuratObject(counts = assay(inSCE, useAssay))
    if ("seuratNormalizedData" %in% names(assays(inSCE))) {
        seuratObject@assays$RNA@data <- assay(inSCE, "seuratNormalizedData")
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if ("seuratScaledData" %in% names(assays(inSCE))) {
        if (!inherits(assay(inSCE, "seuratScaledData"), "matrix")) {
            seuratObject@assays$RNA@scale.data <- as.matrix(assay(inSCE, "seuratScaledData"))
        }
        else {
            seuratObject@assays$RNA@scale.data <- assay(inSCE, "seuratScaledData")
        }
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && length(inSCE@metadata[["seurat"]]@assays$RNA@var.features) > 0) {
        seuratObject@assays$RNA@var.features <- inSCE@metadata[["seurat"]]@assays$RNA@var.features
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@reductions$pca)) {
        seuratObject@reductions$pca <- inSCE@metadata[["seurat"]]@reductions$pca
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@assays$RNA@meta.features)) {
        seuratObject@assays$RNA@meta.features <- inSCE@metadata[["seurat"]]@assays$RNA@meta.features
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@reductions$ica)) {
        seuratObject@reductions$ica <- inSCE@metadata[["seurat"]]@reductions$ica
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@reductions$tsne)) {
        seuratObject@reductions$tsne <- inSCE@metadata[["seurat"]]@reductions$tsne
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@reductions$umap)) {
        seuratObject@reductions$umap <- inSCE@metadata[["seurat"]]@reductions$umap
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@meta.data)) {
        seuratObject@meta.data <- inSCE@metadata[["seurat"]]@meta.data
    }
    if (!is.null(inSCE@metadata[["seurat"]]) && !is.null(inSCE@metadata[["seurat"]]@commands)) {
        seuratObject@commands <- inSCE@metadata[["seurat"]]@commands
    }
    return(seuratObject)
}

#' seuratFindClusters
#' Computes the clusters from the input sce object and stores them back in sce object
#' @param inSCE (sce) object from which clusters should be computed and stored in
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting#' 
#' @param reduction selected reduction method to use for computing clusters ("pca" or "ica", default is "pca")
#' @param dims numeric value of how many components to use for computing clusters (default is 10)
#' @param algorithm selected algorithm to compute clusters (default is "original Louvain algorithm")
#' @param groupSingletons boolean if singletons should be grouped together or not (TRUE or FALSE, default is TRUE)
#' @return sceObject; updated sce object which now contains the computed clusters
#' @export
seuratFindClusters <- function(inSCE, useAssay, geneNamesSeurat, reduction, dims, algorithm, groupSingletons) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    seuratObject <- Seurat::FindNeighbors(seuratObject, reduction = reduction, dims = seq(dims))
    no_algorithm <- 1
    if (algorithm == "original Louvain algorithm") {
        no_algorithm = 1
    } else if (algorithm == "Louvain algorithm with multilevel refinement") {
        no_algorithm = 2
    } else if (algorithm == "SLM algorithm") {
        no_algorithm = 3
    }
    seuratObject <- Seurat::FindClusters(seuratObject, algorithm = no_algorithm, group.singletons = groupSingletons)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratRunTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back into the sce object
#' @param inSCE (sce) object on which to compute the tSNE
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction selected reduction algorithm to use for computing tSNE ("pca" or "ica", default is "pca")
#' @param dims numerical value of how many reduction components to use for tSNE computation (default is 10)
#' @return Updated sce object with tSNE computations stored
#' @export
seuratRunTSNE <- function(inSCE, useAssay, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    seuratObject <- Seurat::RunTSNE(seuratObject, reduction = reduction, dims = 1:dims)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' seuratRunUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back into the sce object
#' @param inSCE (sce) object on which to compute the UMAP
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction selected reduction algorithm to use for computing UMAP ("pca" or "ica", default is "pca")
#' @param dims numerical value of how many reduction components to use for UMAP computation (default is 10)
#' @return Updated sce object with UMAP computations stored
#' @export
seuratRunUMAP <- function(inSCE, useAssay, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    seuratObject <- Seurat::RunUMAP(seuratObject, reduction = reduction, dims = 1:dims)
    inSCE <- .addSeuratToMetaDataSCE(inSCE, seuratObject)
    return(inSCE)
}

#' .seuratGetVariableFeatures
#' Retrieves the requested number of variable feature names
#' @param inSCE (sce) object from which to extract the variable feature names
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param numberOfFeatures numerical value indicating how many feature names should be retrieved (default is 100)
#' @return list() of variable feature names
.seuratGetVariableFeatures <- function(inSCE, useAssay, geneNamesSeurat, numberOfFeatures) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    if (length(seuratObject@assays$RNA@var.features) > 0) {
        return(seuratObject@assays$RNA@var.features[1:numberOfFeatures])
    }
}

#' seuratElbowPlot
#' Computes the plot object for elbow plot from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute the elbow plot (pca should be computed)
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param significantPC a numerical value indicating the number of significant principal components (used to alter the color of the significant components)
#' @return plot object
#' @export
seuratElbowPlot <- function(inSCE, useAssay, geneNamesSeurat, significantPC) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    plot <- Seurat::ElbowPlot(seuratObject)
    plot <- ggplot2::ggplot_build(plot)
    for (i in seq(significantPC)) {
        plot$data[[1]]$shape[i] <- 16
        plot$data[[1]]$colour[i] <- "red"
        plot$data[[1]]$size[i] <- 3.5
    }
    plot <- ggplot2::ggplot_gtable(plot)
    plot <- ggplotify::as.ggplot(plot)
    return(plot)
}

#' seuratJackStrawPlot
#' Computes the plot object for jackstraw plot from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute the jackstraw plot (pca should be computed)
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims a numerical value indicating how many components to use in the computation of jackstraw plot from pca (default is number of pca components computed)
#' @return plot object
#' @export
seuratJackStrawPlot <- function(inSCE, useAssay, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    return(Seurat::JackStrawPlot(seuratObject, dims = 1:dims))
}

#' seuratComputeHeatmap
#' Computes the heatmap plot object from the pca slot in the input sce object
#' @param inSCE (sce) object from which to compute heatmap (pca should be computed)
#' @param useAssay which assay to use from sce object
#' @param geneNamesSeurat a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims numerical value indicating how many components to use for the computation of heatmap plot object (default is number of pca components computed) 
#' @return plot object
#' @export
seuratComputeHeatmap <- function(inSCE, useAssay, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat)
    return(Seurat::DimHeatmap(seuratObject, dims = 1:dims, fast = FALSE, combine = FALSE))
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


#' seuratSCTransform
#' Runs the SCTransform function from seurat library to transform the input data
#' @param inSCE the input sce object containing the counts
#' @param newAssayName name for the transformed counts assay (default is "SCTCounts")
#' @param useAssay which counts assay to use from input sce object for transformation
#' @param geneNamesSeurat genenames for consistent formatting
#' @export
seuratSCTransform <- function(inSCE, newAssayName = "SCTCounts", useAssay, geneNamesSeurat) {
    seuratObject <- Seurat::SCTransform(
        object = convertSCEToSeurat(inSCE, useAssay, geneNamesSeurat),
        assay = "RNA",
        new.assay.name = "SCTransform",
        do.correct.umi = FALSE)
    inSCE <- .updateAssaySCE(inSCE, geneNamesSeurat, seuratObject, newAssayName, "data", slotSeurat = "SCTransform")
    return(inSCE)
}

# ----