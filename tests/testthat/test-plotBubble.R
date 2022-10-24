library(singleCellTK)
library(ggplot2)

test_that(desc = "Testing plotBubble.R", {
  data("pbmc3k_2.7.1_sce")
  bubblePlot <- plotBubble(pbmc3k_2.7.1_sce, gene="CD3E", title="cluster test")
  metrics <- runClusterSummaryMetrics(pbmc3k_2.7.1_sce, useAssay="logcounts", "IL7R", clusters="cluster")
  ggBubblePlot <- .ggBubble(metrics)
  
  testthat::expect_true(is.ggplot(bubblePlot))
  testthat::expect_true(is.ggplot(ggBubblePlot))
})