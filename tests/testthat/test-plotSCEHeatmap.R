library(singleCellTK)
context("Testing plotSCEHeatmap")

data("scExample", package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
sce <- scaterlogNormCounts(sce, "logcounts")
rowData(sce)$fgrp <- rep(c("f1", "f2", "f3"), length.out = nrow(sce))
features <- rownames(sce)[1:30]

heatmapMatrix <- function(hm) {
  if (methods::is(hm, "HeatmapList")) hm <- hm@ht_list[[1]]
  hm@matrix
}

test_that("aggregateCol averages cells per group without deprecations", {
  expect_no_warning(
    hm <- plotSCEHeatmap(sce, useAssay = "logcounts",
                         featureIndex = features, aggregateCol = "type",
                         scale = FALSE, trim = NULL),
    class = "deprecatedWarning"
  )
  mat <- heatmapMatrix(hm)
  logc <- as.matrix(assay(sce, "logcounts")[features, ])
  expected <- vapply(sort(unique(sce$type)),
                     function(g) rowMeans(logc[, sce$type == g, drop = FALSE]),
                     numeric(length(features)))
  expect_equal(colnames(mat), colnames(expected))
  expect_equal(unname(mat), unname(expected))
})

test_that("aggregateRow averages features per group", {
  hm <- plotSCEHeatmap(sce, useAssay = "logcounts",
                       featureIndex = features, aggregateRow = "fgrp",
                       scale = FALSE, trim = NULL)
  mat <- heatmapMatrix(hm)
  logc <- as.matrix(assay(sce, "logcounts")[features, ])
  grp <- rowData(sce)[features, "fgrp"]
  expected <- rowsum(logc, grp) / as.vector(table(grp))
  expect_equal(unname(mat), unname(expected))
})
