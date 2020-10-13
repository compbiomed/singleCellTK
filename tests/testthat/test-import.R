library('SummarizedExperiment')
library('SingleCellExperiment')
library('Seurat')
library('testthat')
library('singleCellTK')
library('GSEABase')

context("Testing import functions")

#####################################
## Importing scRNA-seq Data Functions
#####################################

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
  if (!reticulate::py_module_available("anndata")){
    skip("'anndata' not available. Skipping testing importAnnData")
  }
  sce <- importAnnData(sampleDirs = system.file("extdata/annData_pbmc_3k", package = "singleCellTK"),
                       sampleNames = 'pbmc3k_20by20')
  expect_true(validObject(sce))
})
##################################
## Importing Gene Set Functions
##################################

test_that(desc = "Testing importGeneSetFromGMT", {
  data(scExample)
  gmt <- system.file("extdata/mito_subset.gmt", package = "singleCellTK")
  sce <- importGeneSetsFromGMT(inSCE = sce, file = gmt, by = NULL,
                               collectionName = "test")
  expect_true(inherits(sce@metadata$sctk$genesets$test[[1]], "GeneSet"))
}) 

test_that(desc = "Testing importGeneSetFromList", {
  data(scExample)
  gs1 <- rownames(sce)[1:10]
  gs2 <- rownames(sce)[11:20]
  gs <- list("geneset1" = gs1, "geneset2" = gs2)
  sce <- importGeneSetsFromList(inSCE = sce,
                                geneSetList = gs,
                                collectionName = "test",
                                by = "rownames")
  expect_true(inherits(sce@metadata$sctk$genesets$test[[1]], "GeneSet"))
}) 

test_that(desc = "Testing importGeneSetFromCollection", {
  data(scExample)
  gs1 <- GSEABase::GeneSet(setName = "geneset1", geneIds = rownames(sce)[1:10])
  gs2 <- GSEABase::GeneSet(setName = "geneset2", geneIds = rownames(sce)[11:20])
  gsc <- GSEABase::GeneSetCollection(list(gs1, gs2))
  sce <- importGeneSetsFromCollection(inSCE = sce,
                                      geneSetCollection = gsc,
                                      collectionName = "test",
                                      by = "rownames")
  expect_true(inherits(sce@metadata$sctk$genesets$test[[1]], "GeneSet"))
}) 

test_that(desc = "Testing importGeneSetFromMSigDB", {
  data(scExample)
  sce <- importGeneSetsFromMSigDB(inSCE = sce,
                                  categoryIDs = "C2-CP",
                                  species = "Homo sapiens",
                                  mapping = "gene_symbol",
                                  by = "feature_name")
  expect_true(inherits(sce@metadata$sctk$genesets$"C2-CP"[[1]], "GeneSet"))
}) 
