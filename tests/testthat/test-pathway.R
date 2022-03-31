# Pathway Analysis
library(singleCellTK)
context("Testing pathway analysis functions")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
sce <- scaterlogNormCounts(sce, assayName = "logcounts")
gs1 <- rownames(sce)[seq(10)]
gs2 <- rownames(sce)[seq(11,20)]
gs <- list("geneset1" = gs1, "geneset2" = gs2)


sce <- importGeneSetsFromList(inSCE = sce, geneSetList = gs, by = "rownames")
sce <- runVAM(inSCE = sce, geneSetCollectionName = "GeneSetCollection", 
              useAssay = "logcounts")

test_that(desc = "Testing import genesets", {
    sce <- importGeneSetsFromList(inSCE = sce, 
                                  geneSetList = gs, 
                                  by = "rownames")
    testthat::expect_equal(sctkListGeneSetCollections(sce), "GeneSetCollection")
})

test_that(desc = "Testing VAM", {
    sce <- runVAM(inSCE = sce, geneSetCollectionName = "GeneSetCollection", 
                  useAssay = "logcounts")
    testthat::expect_true("VAM_GeneSetCollection_Distance" %in% reducedDimNames(sce))
    testthat::expect_true("VAM_GeneSetCollection_CDF" %in% reducedDimNames(sce))
    
    vamPlot <- plotPathway(sce, resultName = "VAM_GeneSetCollection_CDF", 
                           geneset = "geneset1", groupBy = "type", 
                           summary = "mean", title = "vam")
    testthat::expect_is(vamPlot, "ggplot")
})

test_that(desc = "Testing GSVA", {
    sce <- runGSVA(sce, geneSetCollectionName = "GeneSetCollection")
    testthat::expect_true("GSVA_GeneSetCollection_Scores" %in% reducedDimNames(sce))
    
    gsvaPlot <- plotPathway(sce, resultName = "GSVA_GeneSetCollection_Scores", 
                            geneset = "geneset1", groupBy = "type", 
                            boxplot = TRUE, violin = FALSE, dots = FALSE)
    testthat::expect_is(gsvaPlot, "ggplot")
})
