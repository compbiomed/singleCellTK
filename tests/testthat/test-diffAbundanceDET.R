# Differential Abundance Tests
library(singleCellTK)
context("Testing differential abundance functions")
data(sceBatches, package = "singleCellTK")

test_that(desc = "Testing diffAbundanceFET", {
    
    sceBatches <- diffAbundanceFET(sceBatches,
                                   cluster = "batch", variable = "cell_type", 
                                   control = "alpha", case = "beta", 
                                   analysisName = "daTest")
    result <- getDiffAbundanceResults(sceBatches, "daTest")
    
    testthat::expect_is(result, "data.frame")
    testthat::expect_equal(nrow(result), length(unique(sceBatches$batch)))
    testthat::expect_equal(ncol(result), 12)
})

test_that(desc = "Testing diffAbundanceFET", {
    plot1 <- plotClusterAbundance(sceBatches,
                                  cluster = "batch", 
                                  variable = "cell_type",
                                  combinePlot = "all")
    testthat::expect_is(plot1, "ggplot")
    
    plot2 <- plotClusterAbundance(sceBatches,
                                  cluster = "batch", 
                                  variable = "cell_type",
                                  combinePlot = "none")
    testthat::expect_is(plot2, "list")
    testthat::expect_is(plot2[[1]], "ggplot")
})
