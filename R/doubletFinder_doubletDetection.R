.runDoubletFinder <- function(counts, seurat.pcs = 30, seurat.res = 2, 
                              seurat.nfeatures = 2000, verbose = FALSE){

        ## Convert to sparse matrix if not already in that format
        counts <- methods::as(counts, "dgCMatrix")
        
        seurat <- Seurat::CreateSeuratObject(counts = counts,
            project = "seurat", min.features = 0)
        seurat <- Seurat::NormalizeData(object = seurat,
            normalization.method="LogNormalize", scale.factor = 10000,
            verbose = verbose)

        seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst",
            nfeatures = seurat.nfeatures, verbose = verbose)

        allGenes <- rownames(seurat)
        seurat <- Seurat::ScaleData(seurat, features = allGenes, verbose = verbose)

        numPc <- min(ncol(seurat@assays$RNA@scale.data) - 1, 50)
        seurat <- Seurat::RunPCA(seurat, features =
                Seurat::VariableFeatures(object = seurat),
                npcs = numPc, verbose = verbose)

        seurat <- Seurat::FindNeighbors(seurat, dims = seurat.pcs, verbose = verbose)
        seurat <- Seurat::FindClusters(seurat, resolution = seurat.res, verbose = verbose)

        sweepResListSeurat <- DoubletFinder::paramSweep_v3(seurat,
            PCs = seurat.pcs, sct = FALSE)
        sweepStatsSeurat <- DoubletFinder::summarizeSweep(sweepResListSeurat,
            GT = FALSE)
        bcmvnSeurat <- DoubletFinder::find.pK(sweepStatsSeurat)
        pkOptimal <- as.numeric(as.matrix(bcmvnSeurat$pK[
            which.max(bcmvnSeurat$MeanBC)]))

        annotations <- seurat@meta.data$seurat_clusters
        homotypicProp <- DoubletFinder::modelHomotypic(annotations)
        nExpPoi <- round(0.075*ncol(seurat@assays$RNA))
        seurat <- DoubletFinder::doubletFinder_v3(seurat,
            PCs = seurat.pcs, pN = 0.25, pK = pkOptimal, nExp = nExpPoi,
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
#' @param seurat.nfeatures Integer. Number of highly variable genes to use. Default 2000.
#' @param seurat.pcs Numeric vector. The PCs used in seurat function to 
#'   determine number of clusters. Default 1:15.
#' @param seurat.res Numeric vector. The resolution parameter used in seurat 
#'   which adjusts the number of clusters determined via the algorithm. Default 2.
#' @param verbose Boolean. Wheter to print messages from Seurat and DoubletFinder. Default FALSE.
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runDoubletFinder(sce_chcl)
#' @export
runDoubletFinder <- function(sce, sample = NULL, seed = 12345, seurat.pcs = 1:15,
                             seurat.res = 2, seurat.nfeatures = 2000, verbose = FALSE){

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

        result <- withr::with_seed(seed, 
                  .runDoubletFinder(counts = sceCountsSub, 
			             seurat.pcs = seurat.pcs,
			             seurat.res = seurat.res,
			             seurat.nfeatures = seurat.nfeatures,
			             verbose = verbose))
    
     doubletScore[sceSubIx] <- result@meta.data$doubletFinderAnnScore
     doubletLabel[sceSubIx] <- result@meta.data$doubletFinderLabel
    }
    colData(sce)$doubletfinder_doublet_score <- doubletScore
    colData(sce)$doubletfinder_doublet_label <- doubletLabel
    return(sce)
}
