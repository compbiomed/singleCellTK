# scran_modelGeneVar.R
library(singleCellTK)
context("Testing scran_modelGeneVar.R")
data(scExample, package = "singleCellTK")
sce <- runSeuratNormalizeData(sce)

test_that(desc = "Testing FindHVG", {
    # Running each method
    sce <- runModelGeneVar(sce, "seuratNormData")
    metricNames.mgv <- metadata(sce)$sctk$runFeatureSelection$modelGeneVar$rowData
    testthat::expect_true(all(metricNames.mgv %in% names(rowData(sce))))

    sce <- runFeatureSelection(sce, "counts", "vst")
    metricNames.vst <- metadata(sce)$sctk$runFeatureSelection$vst$rowData
    testthat::expect_true(all(metricNames.vst %in% names(rowData(sce))))

    sce <- runFeatureSelection(sce, "seuratNormData", "dispersion")
    metricNames.disp <- metadata(sce)$sctk$runFeatureSelection$dispersion$rowData
    testthat::expect_true(all(metricNames.disp %in% names(rowData(sce))))

    sce <- runFeatureSelection(sce, "seuratNormData", "mean.var.plot")
    metricNames.mvp <- metadata(sce)$sctk$runFeatureSelection$mean.var.plot$rowData
    testthat::expect_true(all(metricNames.mvp %in% names(rowData(sce))))

    # Test accessor functions
    sce <- setTopHVG(sce, "modelGeneVar", hvgNumber = 2000, altExp = TRUE, featureSubsetName = NULL)
    nHVG <- length(getTopHVG(sce, useFeatureSubset = "HVG_modelGeneVar2000"))
    
    # testthat::expect_true(is.logical(rowData(sce)$HVG_modelGeneVar2000))
    # testthat::expect_equal(nrow(altExp(sce, "HVG_modelGeneVar2000")), nHVG)
    # testthat::expect_equal(metadata(sce)$sctk$featureSubsets$HVG_modelGeneVar2000$useAssay,
    #                        "seuratNormData")
    # 
    # hvgs <- getTopHVG(sce, "mean.var.plot", hvgNumber = 2000,
    #                   featureDisplay = "feature_name")
    # testthat::expect_false(all(startsWith(hvgs, "ENSG00000")))
    # hvgs <- getTopHVG(sce, "vst", hvgNumber = 2000)
    # vm1 <- plotTopHVG(sce, "dispersion")
    # vm2 <- plotTopHVG(sce, "modelGeneVar", hvgNumber = 30, labelsCount = 10,
    #                  featureDisplay = "feature_name")
    # vm3 <- plotTopHVG(sce, method = "mean.var.plot",
    #                   useFeatureSubset = "HVG_modelGeneVar2000")
    # testthat::expect_true(inherits(vm1, "ggplot"))
    # testthat::expect_true(inherits(vm2, "ggplot"))
    # testthat::expect_true(inherits(vm3, "ggplot"))
})
