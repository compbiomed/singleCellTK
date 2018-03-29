context("differential expression tests")

test_that("scDiffEx and scDiffEx_limma() functions should give the same result", {
  #limma
  # two levels, no covariates
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
  # two levels, covariates
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level1class",
                  covariates = "tissue",
                  ntop = nrow(nozeros),
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level1class",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, biomarker
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "biomarker",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "biomarker",
                         levelofinterest = "Oligo6")
  expect_equal(res, res2)
  # nine levels, covariates, biomarker
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "biomarker",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6",
                  covariates = "tissue")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "biomarker",
                         levelofinterest = "Oligo6",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "coef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "coef",
                         levelofinterest = "Oligo6")
  expect_equal(res, res2)
  # nine levels, covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "coef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6",
                  covariates = "tissue")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "coef",
                         levelofinterest = "Oligo6",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "allcoef",
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "allcoef")
  expect_equal(res, res2)
  # nine levels, covariates, allcoef
  res <- scDiffEx(inSCESet = nozeros,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysis_type = "allcoef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  covariates = "tissue")
  res2 <- scDiffEx_limma(inSCESet = nozeros,
                         use_assay = "logcounts",
                         condition = "level2class",
                         analysis_type = "allcoef",
                         covariates = "tissue")
  expect_equal(res, res2)
  expect_error(scDiffEx(inSCESet = nozeros,
                        use_assay = "logcounts",
                        condition = "level2class",
                        ntop = nrow(nozeros),
                        analysis_type = "false",
                        usesig = FALSE,
                        diffexmethod = "limma",
                        covariates = "tissue"),
               "Unrecognized analysis type, false")
  expect_error(scDiffEx_limma(inSCESet = nozeros,
                              use_assay = "logcounts",
                              condition = "level2class",
                              analysis_type = "false",
                              covariates = "tissue"),
               "Unrecognized analysis type, false")
  rm(nozeros)
})

test_that("scDiffEx and scDiffEx_deseq2() functions should give the same result", {
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
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level1class",
                  ntop = nrow(subset),
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level1class", covariates = "tissue")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, no covariates
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysis_type = "biomarker",
                  levelofinterest = "Oligo6",
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class",
                          analysis_type = "biomarker",
                          levelofinterest = "Oligo6")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, covariates
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysis_type = "biomarker",
                  levelofinterest = "Oligo6",
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysis_type = "biomarker",
                          levelofinterest = "Oligo6")
  expect_equal(res, res2[rownames(res), ])
  #contrast, no covariates
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysis_type = "contrast",
                  levelofinterest = "Oligo6",
                  controlLevel = "Oligo5",
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class",
                          analysis_type = "contrast",
                          levelofinterest = "Oligo6",controlLevel = "Oligo5")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, covariates
  res <- scDiffEx(inSCESet = subset,
                  use_assay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysis_type = "contrast",
                  levelofinterest = "Oligo6",
                  controlLevel = "Oligo5",
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysis_type = "contrast",
                          levelofinterest = "Oligo6",controlLevel = "Oligo5")
  expect_equal(res, res2[rownames(res), ])
  #full-reduced, no covariates
  res <- scDiffEx(inSCESet = subset, use_assay = "counts",
                  condition = "level2class", ntop = nrow(subset),
                  analysis_type = "fullreduced", usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class",
                          analysis_type = "fullreduced")
  expect_equal(res, res2[rownames(res), ])
  #full-reduced, covariates
  res <- scDiffEx(inSCESet = subset, use_assay = "counts",
                  condition = "level2class", ntop = nrow(subset),
                  analysis_type = "fullreduced", usesig = FALSE,
                  diffexmethod = "DESeq2", covariates = "tissue")
  res2 <- scDiffEx_deseq2(inSCESet = subset, use_assay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysis_type = "fullreduced")
  expect_equal(res, res2[rownames(res), ])
  rm(subset)
})

test_that("scDiffEx and scDiffEx_anova() functions should give the same result", {
  #anova
  res <- scDiffEx(inSCESet = mouse_brain_subset_sce,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(mouse_brain_subset_sce),
                  usesig = FALSE,
                  diffexmethod = "ANOVA")
  res2 <- scDiffEx_anova(mouse_brain_subset_sce, use_assay = "logcounts", condition = "level2class")
  expect_equal(res, res2[rownames(res), ])

  #anova covariates
  res <- scDiffEx(inSCESet = mouse_brain_subset_sce,
                  use_assay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(mouse_brain_subset_sce),
                  usesig = FALSE,
                  diffexmethod = "ANOVA",
                  covariates = "tissue")
  res2 <- scDiffEx_anova(mouse_brain_subset_sce, use_assay = "logcounts",
                         condition = "level2class", covariates = "tissue")
  expect_equal(res, res2[rownames(res), ])
})

test_that("condition and covariates should accept numeric vectors", {
  #deseq2
  subset <- mouse_brain_subset_sce[rownames(mouse_brain_subset_sce)[order(rowSums(assay(mouse_brain_subset_sce, "counts")), decreasing = TRUE)][1:100], ]
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                          condition = "age", ntop = nrow(subset),
                          usesig = FALSE, diffexmethod = "DESeq2"), NA)
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                        condition = "level1class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "DESeq2"), NA)
  #limma
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                        condition = "age", ntop = nrow(subset),
                        usesig = FALSE, diffexmethod = "limma"), NA)
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                        condition = "level1class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "limma"), NA)
  #anova
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                        condition = "age", ntop = nrow(subset),
                        usesig = FALSE, diffexmethod = "ANOVA"), NA)
  expect_error(scDiffEx(inSCESet = subset, use_assay = "counts",
                        condition = "level2class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "ANOVA"), NA)
})
