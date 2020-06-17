library('SummarizedExperiment')
library('SingleCellExperiment')
library('Seurat')
library('testthat')
library('singleCellTK')

context("Testing import functions")

test_that(desc = "Testing importBUStools", {
  
  sce <- importBUStools(BUStoolsDirs = system.file("extdata/BUStools_PBMC_1k_v3_20x20/genecount/", package = "singleCellTK"),
                        samples = "PBMC_1k_v3_20x20")
  expect_true(validObject(sce))
})


test_that(desc = "Testing importCellRanger", {
  sce <- importCellRanger(cellRangerDirs = system.file("extdata",package = "singleCellTK"),
                          sampleDirs = "hgmm_1k_v3_20x20",
                          sampleNames = "hgmm1kv3",
                          dataType = "filtered")
  expect_true(validObject(sce))
})

test_that(desc = "Testing importDropEst", {
  sce <- importDropEst(sampleDirs = system.file("extdata/dropEst_scg71",package = "singleCellTK"),
                       sampleNames = 'scg71')
  expect_true(validObject(sce))
})


test_that(desc = "Testing importSeqc", {
  sce <- importSEQC(seqcDirs = system.file("extdata/pbmc_1k_50x50",package = "singleCellTK"),
                    samples = "pbmc_1k_50x50",
                    prefix = "pbmc_1k",
                    combinedSample = FALSE)
  expect_true(validObject(sce))
})

test_that(desc = "Testing importSTARSolo", {
  sce <- importSTARsolo(STARsoloDirs = system.file("extdata/STARsolo_PBMC_1k_v3_20x20",package = "singleCellTK"),
                        samples = "PBMC_1k_v3_20x20")
  
  expect_true(validObject(sce))
})

test_that(desc = "Testing importOptimus", {
  if (!reticulate::py_module_available("scipy.sparse") || (!reticulate::py_module_available("numpy"))){
    skip("scipy.sparse or numpy not available. Skipping testing importOptimus")
  }
  sce <- importOptimus(OptimusDirs = system.file("extdata/Optimus_20x1000",package = "singleCellTK"),
                       samples = "Optimus_20x1000")
  expect_true(validObject(sce))
}) 

test_that(desc = "Testing importAnnData", {
  sce <- importAnnData(sampleDirs = system.file("extdata/annData_pbmc_3k", package = "singleCellTK"),
                       sampleNames = 'pbmc3k_20by20')
  expect_true(validObject(sce))
})