library(singleCellTK)
sce <- importExampleData("pbmc3k")

# Run QC
sce <- runCellQC(sce, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                 mitoRef = "human", mitoIDType = "symbol", 
                 mitoGeneLocation = "rownames", seed = 12345)

sce <- subsetSCECols(sce, colData = c("total > 600", 
                                      "detected > 300",
                                      "mito_percent < 5"))

sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                        normalizationMethod = "logNormCounts")

sce <- runModelGeneVar(sce, useAssay = "logcounts")
sce <- setTopHVG(sce, method = "modelGeneVar", hvgNumber = 2000, featureSubsetName = "hvg2000")
hvg <- getTopHVG(sce, useFeatureSubset = "hvg2000")

sce <- scaterPCA(sce, useFeatureSubset = "hvg2000", seed = 12345)
sce <- runUMAP(sce, useReducedDim = "PCA", initialDims = 10, seed = 12345)
sce <- runScranSNN(sce, "PCA", clusterName = "cluster", nComp = 10,
                   weightType = "jaccard", k = 14, seed = 12345)

sce <- runFindMarker(sce, useAssay = "logcounts", method = "wilcox", cluster = "cluster")


topMarkers <- getFindMarkerTopTable(sce, topN = 1, log2fcThreshold = 0, fdrThreshold = 0.05,
                                    minClustExprPerc = 0.5, maxCtrlExprPerc = 0.5,
                                    minMeanExpr = 0)

pbmc3k_2.7.1_sce <- runSingleR(sce, useAssay = "logcounts", level = "fine")
usethis::use_data(pbmc3k_2.7.1_sce, overwrite = TRUE)
