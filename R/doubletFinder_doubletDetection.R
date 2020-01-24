.runDoubletFinder <- function(counts, ...){
        seurat <- Seurat::CreateSeuratObject(counts = counts,
            project = "seurat", min.features = 0)
        seurat <- Seurat::NormalizeData(object = seurat,
            normalization.method="LogNormalize", scale.factor = 10000,
            verbose = F)

        seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst",
            nfeatures = 2000, verbose = F)

        allGenes <- rownames(seurat)
        seurat <- Seurat::ScaleData(seurat, features = allGenes, verbose = F)

        numPc <- min(ncol(seurat@assays$RNA@scale.data) - 1, 50)
        seurat <- Seurat::RunPCA(seurat, features =
                Seurat::VariableFeatures(object = seurat),
                npcs = numPc, verbose = F)

        seurat <- Seurat::FindNeighbors(seurat, dims = 1:15, verbose = F)
        seurat <- Seurat::FindClusters(seurat, resolution = 1.2, verbose = F)

        sweepResListSeurat <- DoubletFinder::paramSweep_v3(seurat,
            PCs = 1:15, sct = FALSE)
        sweepStatsSeurat <- DoubletFinder::summarizeSweep(sweepResListSeurat,
            GT = FALSE)
        bcmvnSeurat <- DoubletFinder::find.pK(sweepStatsSeurat)
        pkOptimal <- as.numeric(as.matrix(bcmvnSeurat$pK[
            which.max(bcmvnSeurat$MeanBC)]))

        annotations <- seurat@meta.data$seurat_clusters
        homotypicProp <- DoubletFinder::modelHomotypic(annotations)
        nExpPoi <- round(0.075*ncol(seurat@assays$RNA))
        seurat <- DoubletFinder::doubletFinder_v3(seurat,
            PCs = 1:15, pN = 0.25, pK = pkOptimal, nExp = nExpPoi,
            reuse.pANN = FALSE, sct = FALSE)

        names(seurat@meta.data)[6] <- "doubletFinderAnnScore"
        names(seurat@meta.data)[7] <- "doubletFinderLabel"
        return(seurat)
}

  


#' @title Generates a doublet score for each cell via doubletFinder
#' @description Uses doubletFinder to determine cells within the dataset
#'  suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @param sample Numeric vector. Each cell will be assigned a sample number.
#' @param seed Seed for the random number generator. Default 12345.
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' sce <- runDoubletFinder(sce)
#' @export
#' @import SummarizedExperiment
#' @import Seurat
#' @import DoubletFinder
runDoubletFinder <- function(sce, sample, seed = 12345, ...){

  if(!is.null(sample)) {
    if(length(sample) != ncol(sce)) {
      stop("'sample' must be the same length as the number of columns in 'sce'")
    }
  } else {
    sample = rep(1, ncol(sce))
  }

  message(paste0(date(), " ... Running 'doubletFinder'"))

    doubletScore <- rep(NA, ncol(sce))
    doubletLabel <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sample))
    
    for(num in allSampleNumbers){
        sceSubIx <- which(sample == num)
        sceCounts <- SingleCellExperiment::counts(sce)
        sceCountsSub <- sceCounts[,sceSubIx]

        result <- suppressMessages(withr::with_seed(seed, 
                  .runDoubletFinder(counts = sceCountsSub, ...)))
    
     doubletScore[sceSubIx] <- result@meta.data$doubletFinderAnnScore
     doubletLabel[sceSubIx] <- result@meta.data$doubletFinderLabel
    }
    colData(sce)$doublet_finder_doublet_score <- doubletScore
    colData(sce)$doublet_finder_doublet_label <- doubletLabel
    return(sce)
}
