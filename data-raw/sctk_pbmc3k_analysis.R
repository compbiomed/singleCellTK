library(singleCellTK)

#################################### README ####################################
# 1. This script has the exactly same contents as 
# `vignette/articles/console_analysis_tutorial.Rmd`, but the result fetching 
# commands (get table, make report, make plots, etc.) are commented.
# 2. It is IMPORTANT to run the `set.seed` command at the right position in 
# order to reproduce the same result as the Tutorial. 

# Import ####
sce <- importExampleData("pbmc3k")
sce <- importMitoGeneSet(sce, reference = "human", id = "symbol",
                         by = "rownames", collectionName = "mito")

# Run QC ####
sce <- runCellQC(sce, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                 collectionName = "mito",
                 geneSetListLocation = "rownames")

# Plotting for QC. Skipped
# plotRunPerCellQCResults(sce)
# plotScDblFinderResults(sce, reducedDimName = "decontX_UMAP")
# plotDecontXResults(sce, reducedDimName = "decontX_UMAP")

sce <- subsetSCECols(sce, colData = c("total > 600", "detected > 300",
                                      "scDblFinder_doublet_score < 0.8",
                                      "subsets_mito_percent < 5"))

# Normalization and Scaling####
sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                        normalizationMethod = "logNormCounts")
sce <- runNormalization(sce, useAssay = "logcounts",
                        outAssayName = "logcounts_scaled", scale = TRUE)

# Variable Feature ####
sce <- runSeuratFindHVG(sce, useAssay = "counts", hvgMethod = "vst")
sce <- getTopHVG(sce, method = "vst", n = 2000, altExp = "hvg")

# Dimensionality Reduction ####
set.seed(12345)
sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
                    useAssay = "logcounts_scaled", useAltExp = "hvg",
                    reducedDimName = "PCA")

# Clustering ####
sce <- runScranSNN(sce, useReducedDim = "PCA", clusterName = "cluster", 
                   algorithm = "louvain", k = 4)

# Embedding ####
sce <- runDimReduce(inSCE = sce, method = "scaterUMAP", 
                    useReducedDim = "PCA", reducedDimName = "UMAP")
# plotSCEDimReduceColData(sce, colorBy = "cluster", reducedDimName = "UMAP")

# findMarker
sce <- findMarkerDiffExp(sce, useAssay = "logcounts", method = "wilcox",
                         cluster = "cluster",
                         log2fcThreshold = 0, fdrThreshold = 0.05,
                         minClustExprPerc = 0, maxCtrlExprPerc = 1,
                         minMeanExpr = 0)
# topMarkers <- findMarkerTopTable(sce, topN = 1, log2fcThreshold = 0, 
#                                  fdrThreshold = 0.05,
#                                  minClustExprPerc = 0.5, 
#                                  maxCtrlExprPerc = 0.4,
#                                  minMeanExpr = 0)
# head(topMarkers, 10)
# markerHm <- plotMarkerDiffExp(sce, topN = 5, log2fcThreshold = 0, 
#                               fdrThreshold = 0.05, minClustExprPerc = 0.5, 
#                               maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
#                               rowLabel = TRUE)

# Differential Expression Analysis ####
sce <- runDEAnalysis(inSCE = sce, method = "wilcox", 
                     class = "cluster", classGroup1 = c(1), classGroup2 = c(5),
                     groupName1 = "cluster1", groupName2 = "cluster5", 
                     analysisName = "cluster1_VS_5")
# deHm <- plotDEGHeatmap(sce, useResult = "cluster1_VS_5", rowLabel = TRUE)

# Cell Type Labeling ####
sce <- runSingleR(sce, useAssay = "logcounts", level = "fine")
plotSCEDimReduceColData(sce, colorBy = "SingleR_hpca_fine_pruned.labels", 
                        reducedDimName = "UMAP", legendSize = 8)

# Pathway Analysis (VAM) ####
sce <- importGeneSetsFromMSigDB(sce, categoryIDs = "H", 
                                species = "Homo sapiens")
sce <- runVAM(sce, geneSetCollectionName = "H", useAssay = "logcounts")
# hallmark <- "HALLMARK_INFLAMMATORY_RESPONSE"
# plotSCEViolin(sce, slotName = "reducedDims", itemName = "VAM_H_CDF", 
#               dimension = hallmark, boxplot = FALSE,
#               xlab = "cluster", groupBy = sce$cluster, 
#               ylab = hallmark)

saveRDS(sce, "sctk_pbmc3k_analysis.rds")
