.runDoubletFinder <- function(counts, seuratPcs, seuratRes, formationRate,
                              seuratNfeatures, verbose = FALSE){

        ## Convert to sparse matrix if not already in that format
        counts <- methods::as(counts, "dgCMatrix")

        seurat <- Seurat::CreateSeuratObject(counts = counts,
            project = "seurat", min.features = 0)
        seurat <- Seurat::NormalizeData(object = seurat,
            normalization.method="LogNormalize", scale.factor = 10000,
            verbose = verbose)

        seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst",
            nfeatures = seuratNfeatures, verbose = verbose)

        allGenes <- rownames(seurat)
        seurat <- Seurat::ScaleData(seurat, features = allGenes, verbose = verbose)

        numPc <- min(ncol(seurat@assays$RNA@scale.data) - 1, 50)
        seurat <- Seurat::RunPCA(seurat, features =
                Seurat::VariableFeatures(object = seurat),
                npcs = numPc, verbose = verbose)

        seurat <- Seurat::FindNeighbors(seurat, dims = seuratPcs, verbose = verbose)
        seurat <- Seurat::FindClusters(seurat, resolution = seuratRes, verbose = verbose)

        invisible(sweepResListSeurat <- DoubletFinder::paramSweep_v3(seurat,
            PCs = seuratPcs, sct = FALSE))
	invisible(sweepStatsSeurat <- DoubletFinder::summarizeSweep(sweepResListSeurat,
            GT = FALSE))
        bcmvnSeurat <- DoubletFinder::find.pK(sweepStatsSeurat)
        pkOptimal <- as.numeric(as.matrix(bcmvnSeurat$pK[
            which.max(bcmvnSeurat$MeanBC)]))
        annotations <- seurat@meta.data$seurat_clusters
        homotypicProp <- DoubletFinder::modelHomotypic(annotations)
        nExpPoi <- round(formationRate*ncol(seurat@assays$RNA))
        seurat <- invisible(DoubletFinder::doubletFinder_v3(seurat,
            PCs = seuratPcs, pN = 0.25, pK = pkOptimal, nExp = nExpPoi,
            reuse.pANN = FALSE, sct = FALSE))
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
#' @param seuratNfeatures Integer. Number of highly variable genes to use. Default 2000.
#' @param assayName  A string specifying which assay in the SCE to use.
#' @param formationRate Doublet formation rate used within algorithm. Default 0.075.
#' @param seuratPcs Numeric vector. The PCs used in seurat function to
#'   determine number of clusters. Default 1:15.
#' @param seuratRes Numeric vector. The resolution parameter used in seurat. Default c(0.5, 1, 1.5, 2).
#'   which adjusts the number of clusters determined via the algorithm. Default 2.
#' @param verbose Boolean. Wheter to print messages from Seurat and DoubletFinder. Default FALSE.
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runDoubletFinder(sce_chcl)
#' @export
runDoubletFinder <- function(sce, sample = NULL, seed = 12345, seuratPcs = 1:15, formationRate = 0.075,
                             seuratRes = c(0.5, 1, 1.5, 2), seuratNfeatures = 2000, assayName = "counts", verbose = FALSE){

  if (!require("DoubletFinder", character.only = TRUE)){
      message("Package `DoubletFinder` not found, cannot be run. Please install 
	via `devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')`")
      return(sce)	
  }


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
    
    
    for (res in seuratRes){
    output <- S4Vectors::DataFrame(row.names = colnames(sce),
            doubletFinder_doublet_score = numeric(ncol(sce)),
	    doubletFinder_doublet_label = numeric(ncol(sce)))

    ## Loop through each sample and run doubletFinder
    samples <- unique(sample)
	
    for (i in seq_len(length(samples))) {
    	sceSampleInd <- sample == samples[i]
    	sceSample <- sce[, sceSampleInd]
	mat <- SummarizedExperiment::assay(sceSample, i = assayName)
	result <- withr::with_seed(seed,
                  .runDoubletFinder(counts = mat,
                                     seuratPcs = seuratPcs,
                                     seuratRes = res,
                                     seuratNfeatures = seuratNfeatures,
                                     formationRate = formationRate,
                                     verbose = verbose))
	output[sceSampleInd, 1] <- result@meta.data$doubletFinderAnnScore
	output[sceSampleInd, 2] <- result@meta.data$doubletFinderLabel
	}
    
    colnames(output) <- paste0(colnames(output), "_Resolution_", res)
    
    colData(sce) = cbind(colData(sce), output)
    }
    return(sce)
}

