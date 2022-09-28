# dimension reduction
library(singleCellTK)
library(testthat)
context("Testing dimension reduction")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
sce <- scaterlogNormCounts(sce, "logcounts")
sce <- runFeatureSelection(sce, useAssay = "counts")
sce <- setTopHVG(sce, featureSubsetName = "hvg")
sce <- setTopHVG(sce, featureSubsetName = "hvgAltExp", altExp = TRUE)

test_that(desc = "Testing scaterPCA", {
    sce <- scaterPCA(sce, useAssay = "logcounts", useFeatureSubset = "hvg",
                     reducedDimName = "PCA1")
    testthat::expect_true("PCA1" %in% reducedDimNames(sce))

    sce <- scaterPCA(sce, useAssay = "hvgAltExplogcounts", useAltExp = "hvgAltExp",
                     reducedDimName = "PCA2")
    testthat::expect_true("PCA2" %in% reducedDimNames(sce))

    expect_error({
        scaterPCA(sce, useAssay = "null", useFeatureSubset = "hvg")
    }, "Specified `useAssay` 'null' not found.")
    expect_error({
        scaterPCA(sce, useAssay = "null", useAltExp = "null", useFeatureSubset = "hvg")
    }, "Specified `useAltExp` 'null' not found.")
    expect_error({
        scaterPCA(sce, useAssay = "null", useAltExp = "hvgAltExp", useFeatureSubset = "hvg")
    }, "Specified `useAssay` 'null' not found in the altExp.")

    p1 <- plotPCA(sce, reducedDimName = "PCA1", colorBy = "type", shape = "type")
    p2 <- plotPCA(sce, reducedDimName = "PCA3", runPCA = TRUE)
    expect_is(p1, "ggplot")
    expect_is(p2, "ggplot")
    expect_error({
        plotPCA(sce, pcX = "foo", reducedDimName = "PCA1")
    }, regexp = "pcX dimension")
    expect_error({
        plotPCA(sce, pcY = "bar", reducedDimName = "PCA1")
    }, regexp = "pcY dimension")
    expect_error({
        plotPCA(sce, reducedDimName = "UMAP")
    }, regexp = "dimension not found")
})

test_that(desc = "Testing scater UMAP", {
    sce <- scaterPCA(sce, useFeatureSubset = "hvg", seed = 12345, reducedDimName = "PCA1")
    sce <- runUMAP(sce, useReducedDim = "PCA1", reducedDimName = "UMAP1")
    testthat::expect_true("UMAP1" %in% reducedDimNames(sce))
    sce <- runUMAP(sce, useAssay = "hvgAltExplogcounts", useReducedDim = NULL,
                   useAltExp = "hvgAltExp",
                   reducedDimName = "UMAP2")
    testthat::expect_true("UMAP2" %in% reducedDimNames(sce))
    # TODO: Still some runable conditions
    expect_error({
        runUMAP(inSCE = 1)
    }, "`inSCE` should be a SingleCellExperiment Object.")
    expect_error({
        runUMAP(sce, useAssay = NULL, useReducedDim = NULL)
    }, "Either `useAssay` or `useReducedDim` has to be specified.")
    expect_error({
        runUMAP(sce, useAltExp = "altexp")
    }, "Specified `useAltExp` 'altexp' not found.")
    expect_error({
        runUMAP(sce, useReducedDim = "TSNE")
    }, "Specified `useReducedDim` 'TSNE' not found.")
    expect_error({
        runUMAP(sce, useAssay = "TSNE", useReducedDim = NULL)
    }, regexp = "Specified `useAssay` 'TSNE' not found.")
    expect_error({
        runUMAP(sce, useReducedDim = "PCA1", sample = "batch")
    }, regexp = "Specified variable 'batch'")
    expect_error({
        runUMAP(sce, useReducedDim = "PCA1", sample = letters)
    }, regexp = "Invalid variable length")

    p1 <- plotUMAP(sce, reducedDimName = "UMAP1", colorBy = "type", shape = "type")
    p2 <- plotUMAP(sce, reducedDimName = "UMAP3", runUMAP = TRUE)
    expect_is(p1, "ggplot")
    expect_is(p2, "ggplot")
    expect_error({
        plotUMAP(sce, reducedDimName = "TSNE")
    }, regexp = "dimension not found")
    reducedDim(sce, "UMAP4") <- cbind(reducedDim(sce, "UMAP1"), reducedDim(sce, "UMAP1"))
    expect_warning({
        plotUMAP(sce, reducedDimName = "UMAP4")
    }, "More than two UMAP dimensions")
})

test_that(desc = "Testing Rtsne TSNE", {
    sce <- scaterPCA(sce, useFeatureSubset = "hvg", seed = 12345, reducedDimName = "PCA1")
    sce <- runTSNE(sce, useReducedDim = "PCA1", reducedDimName = "TSNE1")
    testthat::expect_true("TSNE1" %in% reducedDimNames(sce))
    sce <- runQuickTSNE(sce, useAssay = "hvgAltExpcounts", useAltExp = "hvgAltExp",
                        reducedDimName = "TSNE2", logNorm = TRUE, nTop = 50)
    testthat::expect_true("TSNE2" %in% reducedDimNames(sce))
    # TODO: Still some runable conditions
    expect_error({
        runTSNE(sce, useAssay = NULL, useReducedDim = NULL)
    }, "Either `useAssay` or `useReducedDim` has to be specified.")
    expect_error({
        runTSNE(sce, useAltExp = "altexp")
    }, "Specified `useAltExp` 'altexp' not found.")

    p1 <- plotTSNE(sce, reducedDimName = "TSNE1", colorBy = "type", shape = "type")
    p2 <- plotTSNE(sce, reducedDimName = "TSNE3", runTSNE = TRUE)
    expect_is(p1, "ggplot")
    expect_is(p2, "ggplot")
    expect_error({
        plotTSNE(sce, reducedDimName = "UMAP")
    }, regexp = "dimension not found")
    reducedDim(sce, "TSNE4") <- cbind(reducedDim(sce, "TSNE1"), reducedDim(sce, "TSNE1"))
    expect_warning({
        plotTSNE(sce, reducedDimName = "TSNE4")
    }, "More than two t-SNE dimensions")

    p3 <- plotDimRed(sce, useReduction = "PCA1")
    p4 <- plotDimRed(sce, useReduction = "PCA1", xAxisLabel = "PC_1", yAxisLabel = "PC_2")
    expect_is(p3, "ggplot")
    expect_is(p4, "ggplot")
})
