library('SummarizedExperiment')
library('SingleCellExperiment')
library('Seurat')
library('testthat')
library('singleCellTK')

input_dir <- file.path("./inst/extdata/")
sapply(paste0('./R/',list.files("./R",pattern="*.R")),source)

context("Testing import functions")

test_that(desc = "Testing importBUStools", {
  
  sce <- importBUStools(BUStoolsDirs = file.path(input_dir,"BUStools_PBMC_1k_v3_20x20/genecount/"),
                        samples = "PBMC_1k_v3_20x20")
  expect_true(validObject(sce))
})


test_that(desc = "Testing importCellRanger", {
  sce <- importCellRanger(cellRangerDirs = file.path(input_dir),
                          sampleDirs = "hgmm_1k_v3_20x20",
                          sampleNames = "hgmm1kv3",
                          dataType = "filtered")
  expect_true(validObject(sce))
})

test_that(desc = "Testing importDropEst", {
  sce <- importDropEst(sampleDirs = file.path(input_dir,"dropEst_scg71"),
                       sampleNames = 'scg71')
  expect_true(validObject(sce))
})


test_that(desc = "Testing importSeqc", {
  sce <- importSEQC(seqcDirs = file.path(input_dir,"pbmc_1k_50x50"),
                    samples = "pbmc_1k_50x50",
                    prefix = "pbmc_1k",
                    combinedSample = FALSE)
  expect_true(validObject(sce))
})

test_that(desc = "Testing importSTARSolo", {
  sce <- importSTARsolo(STARsoloDirs = file.path(input_dir,"STARsolo_PBMC_1k_v3_20x20"),
                        samples = "PBMC_1k_v3_20x20")
  
  expect_true(validObject(sce))
})

test_that(desc = "Testing importOptimus", {
  if (!reticulate::py_module_available("scipy.sparse") || (!reticulate::py_module_available("numpy"))){
    skip("scipy.sparse or numpy not available. Skipping testing importOptimus")
  }
  sce <- importOptimus(OptimusDirs = file.path(input_dir,"Optimus_20x1000"),
                       samples = "Optimus_20x1000")
  expect_true(validObject(sce))
}) 