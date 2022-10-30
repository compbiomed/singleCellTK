library(singleCellTK)
library(ggplot2)

test_that(desc = "Testing plotBubble.R", {
  data("scExample")
  bubblePlot <- plotBubble(sce, useAssay="counts", gene="B2M", displayName="feature_name", cluster="type", title="cluster test")
  
  testthat::expect_true(is.ggplot(bubblePlot))
})