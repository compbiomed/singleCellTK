#' A wrapper function run complete Seurat workflow.
#'
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay Specify the assay to use with Seurat workflow.
#'
#' @return \code{SingleCellExperiment} with results from Seurat workflow stored.
#' @export
#'
#' @examples
#' data(sce_chcl, package = "scds")
#' sce_chcl <- runSeurat(inSCE = sce_chcl, useAssay = "counts")
runSeurat <- function(inSCE, useAssay){
  
  message("Normalizing Data")
  # Normalize Data
  inSCE <- seuratNormalizeData(
    inSCE = inSCE,
    useAssay = useAssay
  )
  
  message("Scaling Data")
  # Scale Data
  inSCE <- seuratScaleData(
    inSCE = inSCE
  )
  
  message("Identifying highly variable genes")
  # Feature Selection
  inSCE <- seuratFindHVG(
    inSCE = inSCE
  )
  
  message("Computing reduced dimensions")
  # Dimensionality Reduction
  inSCE <- seuratPCA(
    inSCE = inSCE
  )
  
  message("Computing tSNE/UMAP")
  # tSNE/UMAP
  inSCE <- seuratRunTSNE(
    inSCE = inSCE,
    useReduction = "pca"
  )
  inSCE <- seuratRunUMAP(
    inSCE = inSCE,
    useReduction = "pca"
  )
  
  message("Identifying clusters in data")
  # Clustering
  inSCE <- seuratFindClusters(
    inSCE = inSCE,
    useReduction = "pca",
    algorithm = "louvain"
  )
  
  message("Identifying marker genes in data")
  # Find Markers
  inSCE <- seuratFindMarkers(
    inSCE = inSCE,
    allGroup = "Seurat_louvain_Resolution0.8"
  )
  
  return(inSCE)
}