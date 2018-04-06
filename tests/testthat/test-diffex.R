context("differential expression tests")

test_that("scDiffEx and scDiffExlimma() functions should give the same result", {
  #limma
  # two levels, no covariates
  nozeros <- mouseBrainSubsetSCE[rowSums(assay(mouseBrainSubsetSCE, "counts")) != 0, ]
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level1class",
                  ntop = nrow(nozeros),
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level1class")
  expect_equal(res, res2)
  # two levels, covariates
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level1class",
                  covariates = "tissue",
                  ntop = nrow(nozeros),
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level1class",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, biomarker
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "biomarker",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "biomarker",
                         levelofinterest = "Oligo6")
  expect_equal(res, res2)
  # nine levels, covariates, biomarker
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "biomarker",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6",
                  covariates = "tissue")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "biomarker",
                         levelofinterest = "Oligo6",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "coef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "coef",
                         levelofinterest = "Oligo6")
  expect_equal(res, res2)
  # nine levels, covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "coef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  levelofinterest = "Oligo6",
                  covariates = "tissue")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "coef",
                         levelofinterest = "Oligo6",
                         covariates = "tissue")
  expect_equal(res, res2)
  # nine levels, no covariates, coef
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "allcoef",
                  usesig = FALSE,
                  diffexmethod = "limma")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "allcoef")
  expect_equal(res, res2)
  # nine levels, covariates, allcoef
  res <- scDiffEx(inSCESet = nozeros,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(nozeros),
                  analysisType = "allcoef",
                  usesig = FALSE,
                  diffexmethod = "limma",
                  covariates = "tissue")
  res2 <- scDiffExlimma(inSCESet = nozeros,
                         useAssay = "logcounts",
                         condition = "level2class",
                         analysisType = "allcoef",
                         covariates = "tissue")
  expect_equal(res, res2)
  expect_error(scDiffEx(inSCESet = nozeros,
                        useAssay = "logcounts",
                        condition = "level2class",
                        ntop = nrow(nozeros),
                        analysisType = "false",
                        usesig = FALSE,
                        diffexmethod = "limma",
                        covariates = "tissue"),
               "Unrecognized analysis type, false")
  expect_error(scDiffExlimma(inSCESet = nozeros,
                              useAssay = "logcounts",
                              condition = "level2class",
                              analysisType = "false",
                              covariates = "tissue"),
               "Unrecognized analysis type, false")
  rm(nozeros)
})

test_that("scDiffEx and scDiffExDESeq2() functions should give the same result", {
  #deseq2
  subset <- mouseBrainSubsetSCE[rownames(mouseBrainSubsetSCE)[order(rowSums(assay(mouseBrainSubsetSCE, "counts")), decreasing = TRUE)][1:100], ]
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level1class",
                  ntop = nrow(subset),
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts", condition = "level1class")
  expect_equal(res, res2[rownames(res), ])
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level1class",
                  ntop = nrow(subset),
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level1class", covariates = "tissue")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, no covariates
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysisType = "biomarker",
                  levelofinterest = "Oligo6",
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class",
                          analysisType = "biomarker",
                          levelofinterest = "Oligo6")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, covariates
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysisType = "biomarker",
                  levelofinterest = "Oligo6",
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysisType = "biomarker",
                          levelofinterest = "Oligo6")
  expect_equal(res, res2[rownames(res), ])
  #contrast, no covariates
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysisType = "contrast",
                  levelofinterest = "Oligo6",
                  controlLevel = "Oligo5",
                  usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class",
                          analysisType = "contrast",
                          levelofinterest = "Oligo6",controlLevel = "Oligo5")
  expect_equal(res, res2[rownames(res), ])
  #biomarker, covariates
  res <- scDiffEx(inSCESet = subset,
                  useAssay = "counts",
                  condition = "level2class",
                  ntop = nrow(subset),
                  analysisType = "contrast",
                  levelofinterest = "Oligo6",
                  controlLevel = "Oligo5",
                  usesig = FALSE,
                  diffexmethod = "DESeq2",
                  covariates = "tissue")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysisType = "contrast",
                          levelofinterest = "Oligo6",controlLevel = "Oligo5")
  expect_equal(res, res2[rownames(res), ])
  #full-reduced, no covariates
  res <- scDiffEx(inSCESet = subset, useAssay = "counts",
                  condition = "level2class", ntop = nrow(subset),
                  analysisType = "fullreduced", usesig = FALSE,
                  diffexmethod = "DESeq2")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class",
                          analysisType = "fullreduced")
  expect_equal(res, res2[rownames(res), ])
  #full-reduced, covariates
  res <- scDiffEx(inSCESet = subset, useAssay = "counts",
                  condition = "level2class", ntop = nrow(subset),
                  analysisType = "fullreduced", usesig = FALSE,
                  diffexmethod = "DESeq2", covariates = "tissue")
  res2 <- scDiffExDESeq2(inSCESet = subset, useAssay = "counts",
                          condition = "level2class", covariates = "tissue",
                          analysisType = "fullreduced")
  expect_equal(res, res2[rownames(res), ])
  rm(subset)
})

test_that("scDiffEx and scDiffExANOVA() functions should give the same result", {
  #anova
  res <- scDiffEx(inSCESet = mouseBrainSubsetSCE,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(mouseBrainSubsetSCE),
                  usesig = FALSE,
                  diffexmethod = "ANOVA")
  res2 <- scDiffExANOVA(mouseBrainSubsetSCE, useAssay = "logcounts", condition = "level2class")
  expect_equal(res, res2[rownames(res), ])

  #anova covariates
  res <- scDiffEx(inSCESet = mouseBrainSubsetSCE,
                  useAssay = "logcounts",
                  condition = "level2class",
                  ntop = nrow(mouseBrainSubsetSCE),
                  usesig = FALSE,
                  diffexmethod = "ANOVA",
                  covariates = "tissue")
  res2 <- scDiffExANOVA(mouseBrainSubsetSCE, useAssay = "logcounts",
                        condition = "level2class", covariates = "tissue")
  expect_equal(res, res2[rownames(res), ])
})

test_that("condition and covariates should accept numeric vectors", {
  #deseq2
  subset <- mouseBrainSubsetSCE[rownames(mouseBrainSubsetSCE)[order(rowSums(assay(mouseBrainSubsetSCE, "counts")), decreasing = TRUE)][1:100], ]
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                          condition = "age", ntop = nrow(subset),
                          usesig = FALSE, diffexmethod = "DESeq2"), NA)
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                        condition = "level1class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "DESeq2"), NA)
  #limma
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                        condition = "age", ntop = nrow(subset),
                        usesig = FALSE, diffexmethod = "limma"), NA)
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                        condition = "level1class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "limma"), NA)
  #anova
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                        condition = "age", ntop = nrow(subset),
                        usesig = FALSE, diffexmethod = "ANOVA"), NA)
  expect_error(scDiffEx(inSCESet = subset, useAssay = "counts",
                        condition = "level2class", covariates = "age",
                        ntop = nrow(subset), usesig = FALSE,
                        diffexmethod = "ANOVA"), NA)
})
