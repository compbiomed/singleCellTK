library(singleCellTK)
context("Testing mergeSCColData")
source("C:/Users/Nathan Palamuttam/campbio BU/singleCellTK/R/miscFunctions.R")

data(scExample, package = "singleCellTK")

sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")

colData(sce)$column_name = rownames(colData(sce))
test_that(desc = "Testing renameClusters", {
  renameClusters(inSCE = smallSCE, clusterName = "louvain_0_2", from = c(0), to = c("a"))
  colData(sce)$louvain_0_2
  expect_equal(ncol(colData(sce)) + 1, ncol(colData(mergedsce)))
})



