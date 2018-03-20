context("differential expression tests")

test_that("scDiffEx and method specific scDiffEx functions should give the same result", {
  #limma
  nozeros <- mouse_brain_subset_sce[rowSums(assay(mouse_brain_subset_sce, "counts")) != 0, ]
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level1class",
                  ntop = nrow(nozeros),
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level1class")
  expect_equal(res, res2)
  rm(nozeros)

  #deseq2
  subset <- mouse_brain_subset_sce[rownames(mouse_brain_subset_sce)[order(rowSums(assay(mouse_brain_subset_sce, "counts")), decreasing = TRUE)][1:100], ]
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level1class",
                  ntop = nrow(subset),
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts", condition = "level1class")
  expect_equal(res, res2[rownames(res), ])
  rm(subset)

  #anova
  res <- scDiffEx(inSCESet = mouse_brain_subset_sce,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(mouse_brain_subset_sce),
                  usesig = FALSE,
                  diffexmethod = "ANOVA")
  res2 <- scDiffEx_anova(mouse_brain_subset_sce, use_assay = "logcounts", condition = "level2class")
  expect_equal(res, res2[rownames(res), ])
})
