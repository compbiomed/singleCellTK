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

test_that(desc = "Testing import single-cell data", {
  allImportEntries <- list(
    samples = list()
  )
  allImportEntries$samples$a <- list(
    type = "busTools",
    params = list(BUStoolsDirs = system.file("extdata/BUStools_PBMC_1k_v3_20x20/genecount/", package = "singleCellTK"),
                  samples = "BUStools_20x20")
  )
  allImportEntries$samples$b <- list(
    type = "cellRanger3",
    params = list(cellRangerDirs = system.file("extdata",package = "singleCellTK"),
                  sampleDirs = "hgmm_1k_v3_20x20",
                  sampleNames = "cellRanger3_20x20",
                  dataType = "filtered")
  )
  allImportEntries$samples$c <- list(
    type = "seqc",
    params = list(seqcDirs = system.file("extdata/pbmc_1k_50x50",package = "singleCellTK"),
                  samples = "seqc_50x50",
                  prefix = "pbmc_1k",
                  combinedSample = FALSE)
  )
  allImportEntries$samples$d <- list(
    type = "starSolo",
    params = list(STARsoloDirs = system.file("extdata/STARsolo_PBMC_1k_v3_20x20",package = "singleCellTK"),
                  samples = "STARsolo_20x20")
  )
  allImportEntries$samples$e <- list(
    type = "files",
    params = list(assayFile = system.file("extdata/pbmc_1k_50x50/pbmc_1k_sparse_molecule_counts.mtx",package = "singleCellTK"),
                  annotFile = system.file("extdata/pbmc_1k_50x50/pbmc_1k_sparse_counts_barcodes.csv",package = "singleCellTK"),
                  featureFile = system.file("extdata/pbmc_1k_50x50/pbmc_1k_sparse_counts_genes.csv",package = "singleCellTK"),
                  assayName = "counts")
  )
  expect_warning({
    sce <- importMultipleSources(allImportEntries)
  }, "column names 'cell_barcode.x'")
  expect_true(validObject(sce))
  summary <- as.data.frame(table(sce$sample))
  expect_true(all(summary$Var1 %in% c("BUStools_20x20", "cellRanger3_20x20", "seqc_50x50", "STARsolo_20x20", "sample")))
  expect_length(summary$Freq, 5)
})

test_that(desc = "Testing importDropEst", {
  sce <- importDropEst(sampleDirs = system.file("extdata/dropEst_scg71",package = "singleCellTK"),
                       sampleNames = 'scg71')
  expect_true(validObject(sce))
})

#
# test_that(desc = "Testing importOptimus", {
#   if (!reticulate::py_module_available("scipy.sparse") || (!reticulate::py_module_available("numpy"))){
#     skip("scipy.sparse or numpy not available. Skipping testing importOptimus")
#   }
#   sce <- importOptimus(OptimusDirs = system.file("extdata/Optimus_20x1000",package = "singleCellTK"),
#                        samples = "Optimus_20x1000")
#   expect_true(validObject(sce))
# })

# test_that(desc = "Testing importAnnData", {
#   if (!reticulate::py_module_available("anndata")){
#     skip("'anndata' not available. Skipping testing importAnnData")
#   }
#   sce <- importAnnData(sampleDirs = system.file("extdata/annData_pbmc_3k", package = "singleCellTK"),
#                        sampleNames = 'pbmc3k_20by20')
#   expect_true(validObject(sce))
# })
##################################
## Importing Gene Set Functions
##################################
data(scExample)
test_that(desc = "Testing importGeneSetFromGMT", {
  gmt <- system.file("extdata/mito_subset.gmt", package = "singleCellTK")
  sce <- importGeneSetsFromGMT(inSCE = sce, file = gmt, by = NULL,
                               collectionName = "test")
  expect_true(inherits(sce@metadata$sctk$genesets$test[[1]], "GeneSet"))
})

test_that(desc = "Testing importGeneSetFromList", {
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
  gs1 <- GSEABase::GeneSet(setName = "geneset1", geneIds = rownames(sce)[1:10])
  gs2 <- GSEABase::GeneSet(setName = "geneset2", geneIds = rownames(sce)[11:20])
  gsc <- GSEABase::GeneSetCollection(list(gs1, gs2))
  sce <- importGeneSetsFromCollection(inSCE = sce,
                                      geneSetCollection = gsc,
                                      collectionName = "test",
                                      by = "rownames")
  expect_true(inherits(sce@metadata$sctk$genesets$test[[1]], "GeneSet"))
  gscList <- sctkListGeneSetCollections(sce)
  expect_equal(gscList, "test")
  gsList <- getGenesetNamesFromCollection(sce, "test")
  expect_equal(gsList, c("geneset1", "geneset2"))
})

test_that(desc = "Testing importGeneSetFromMSigDB", {
  sce <- importGeneSetsFromMSigDB(inSCE = sce,
                                  categoryIDs = "C2-CP",
                                  species = "Homo sapiens",
                                  mapping = "gene_symbol",
                                  by = "feature_name")
  expect_true(inherits(sce@metadata$sctk$genesets$"C2-CP"[[1]], "GeneSet"))
})

data(scExample)
test_that(desc = "Testing export", {
  exportSCEtoFlatFile(sce)
  expect_true(file.exists("./assays/SCE_counts.mtx.gz"))
  expect_true(file.exists("./SCE_cellData.txt.gz"))
  expect_true(file.exists("./SCE_featureData.txt.gz"))
})

