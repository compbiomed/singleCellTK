context("check for conditions to visualize plots")

#how to pass an object which is not a sctk object
test_that("Should fail when the object supplied is not a sctk/ sctkexperiment object", {
  expect_error(visPlot(matrix(1:100), "counts", "boxplot", NULL,
                       c("Cmtm5", "C1qa")),
               "Please use a singleCellTK or a SCtkExperiment object")
})

test_that("Should fail when no condition is provided", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "boxplot", NULL,
                       c("Cmtm5", "C1qa")), "Please supply a condition")
})

test_that("Should fail when an incorrect plot name is provided", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "nothingplot",
                       "level1class", c("Cmtm5", "C1qa")),
               "method 'nothingplot' is not a valid method.")
})

test_that("Should fail when an incorrect assay name is provided", {
  expect_error(visPlot(mouseBrainSubsetSCE, "david", "boxplot", "level1class",
                       c("Cmtm5", "C1qa")),
               "assay 'david' does not exist.")
})

test_that("Should fail when an incorrect gene name is provided", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "boxplot", "level1class",
                       "david_gene"),
               "Gene in gene list not found in input object.")
})

test_that("Should fail when more than 1 condition is provided", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "boxplot",
                       c("age", "diameter"), "Cmtm5"),
               "Only 1 condition allowed")
})

test_that("Should fail when the condition is not a factor for boxplot", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "boxplot", "age",
                       "Cmtm5"),
               "Boxplot requires condition to be a factor, use scatterplot instead")
})

test_that("Should fail when there are more than 16 genes for boxplot", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "boxplot", "age",
                       c("Tspan12", "Tshz1", "Fnbp1l",
      "Adamts15", "Cldn12", "Rxfp1", "2310042E22Rik", "Sema3c", "Jam2",
      "Apbb1ip", "Frem2", "BC005764", "Deptor", "C130030K03Rik", "Klhl13",
      "Tnfaip8l3", "Cmtm5")),
      "Maximum limit of genes reached. Please enter 16 or less genes.")
})

test_that("Should fail when the condition is a factor for scatterplot", {
  expect_error(visPlot(mouseBrainSubsetSCE, "logcounts", "scatterplot",
                       "level1class", c("Cmtm5", "C1qa")),
               "Scatterplot requires a condition to be continuous, use boxplot as an alternative")
})

test_that("Should fail when the condition is NULL for scatterplot", {
  expect_error(visPlot(mouseBrainSubsetSCE, "logcounts", "scatterplot", NULL,
                       c("Cmtm5", "C1qa")),
               "Please supply a condition")
})

test_that("Should fail when a condition is supplied for barplot", {
  expect_error(visPlot(mouseBrainSubsetSCE, "logcounts", "barplot", "age",
                       c("Cmtm5", "C1qa")),
               "Barplot doesn't require a condition, use scatterplot or boxplot instead")
})

test_that("Should fail when the input genes supplied for a heatmap has a sum of zero counts", {
  expect_error(visPlot(mouseBrainSubsetSCE, "counts", "heatmap", NULL,
                       "Adamts15"),
               "Gene Adamts15 has zero variance, please filter and continue.")
})
