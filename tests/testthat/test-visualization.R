# dimensionality reduction algorithms
library(singleCellTK)
context("Testing dimensionality reduction algorithms")
data(scExample, package = "singleCellTK")
sceDroplet <- sce
sce <- sce[, colData(sce)$type != 'EmptyDroplet']
sampleVector <- c(rep("Sample1", 100), rep("Sample2", 95))
sceres <- getUMAP(inSCE = sce, useAssay = "counts", logNorm = TRUE, sample = sampleVector, nNeighbors = 10, reducedDimName = "UMAP",
                nIterations = 20, alpha = 1, minDist = 0.01, pca = TRUE, initialDims = 20)

test_that(desc = "Testing getUMAP", {
        expect_equal(names(reducedDims(sceres)), "UMAP")
	expect_equal(nrow(reducedDim(sceres, "UMAP")), ncol(sce))
})

test_that(desc = "Testing plotSCEScatter functions", {
    p1 <- plotSCEScatter(inSCE = sceres, legendTitle = NULL,
        slot = "assays", annotation = "counts", feature = "ENSG00000251562",
        reducedDimName = "UMAP", labelClusters = FALSE, sample = sampleVector)
    expect_is(p1, "ggplot")
    p2 <- plotSCEDimReduceFeatures(inSCE = sceres, feature = "ENSG00000251562",
        shape = NULL, reducedDimName = "UMAP",
        useAssay = "counts", xlab = "UMAP1", ylab = "UMAP2", sample = sampleVector)
    expect_is(p2, "ggplot")
    p3 <- plotSCEDimReduceColData(inSCE = sceres, colorBy = "type",
        shape = NULL, conditionClass = "factor",
        reducedDimName = "UMAP",
        xlab = "UMAP1", ylab = "UMAP2", labelClusters = TRUE, sample = sampleVector)
    expect_is(p3, "ggplot")
})

test_that(desc = "Testing plotSCEViolin functions", {
    p1 <- plotSCEViolin(inSCE = sceres, slot = "assays",
        annotation = "counts", feature = "ENSG00000251562", groupby = "type", sample = sampleVector)
    expect_is(p1, "ggplot")
    p2 <- plotSCEViolinAssayData(inSCE = sceres,
        feature = "ENSG00000251562", groupby = "type", sample = sampleVector)
    expect_is(p2, "ggplot")
    p3 <- plotSCEViolinColData(inSCE = sceres,
        coldata = "type", groupby = "sample", sample = sampleVector)
    expect_is(p3, "ggplot")
})


test_that(desc = "Testing plotResults functions", {
  sceres <- sceres[, colData(sceres)$type != 'EmptyDroplet']
  sceres <- runCellQC(sceres, algorithms = c("QCMetrics", "cxds", "bcds", "cxds_bcds_hybrid",
                                             "scrublet", "doubletFinder", "decontX"))
  sceres <- runDoubletCells(sceres, size.factors.norm = rep(1, ncol(sceres)))

  r1 <- plotRunPerCellQCResults(inSCE = sceres, sample = sampleVector, combinePlot = "all")
    expect_is(r1, c("gg","ggplot"))
  r2 <- plotScrubletResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r2, c("gg","ggplot"))
  r3 <- plotDoubletCellsResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r3, c("gg","ggplot"))
  r4 <- plotDoubletFinderResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r4, c("gg","ggplot"))
  r5 <- plotCxdsResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r5,  c("gg","ggplot"))
  r6 <- plotBcdsResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r6,  c("gg","ggplot"))
  r7 <- plotScdsHybridResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r7,  c("gg","ggplot"))
  r8 <- plotDecontXResults(inSCE = sceres, reducedDimName="UMAP", sample = sampleVector, combinePlot = "all")
    expect_is(r8, c("gg","ggplot"))
  
  sceDroplet <- runDropletQC(sceDroplet)
  r9 <- plotEmptyDropsResults(inSCE = sceDroplet, sample = c(rep("Sample1", 100), rep("Sample2", 290)))
    expect_is(r9, "list")
})

