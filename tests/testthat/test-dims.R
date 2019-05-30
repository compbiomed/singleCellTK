context("reduced dimension functions")

test_that("getPCA", {
  expect_error(getPCA(mouseBrainSubsetSCE, "abcd"),
               "abcd not in the assay list")
  expect_error(getPCA(mouseBrainSubsetSCE[1:100, ], "abcd"),
               "abcd not in the assay list")
  expect_is(getPCA(mouseBrainSubsetSCE[1:100, ], "counts"),
            "SCtkExperiment")
  assay(mouseBrainSubsetSCE, "test") <- data.frame(assay(mouseBrainSubsetSCE,
                                                         "counts"))
  expect_error(getPCA(mouseBrainSubsetSCE[1:100, ], "test"),
               "Input matrix test is not a matrix")
})

test_that("plotPCA", {
  expect_is(plotPCA(mouseBrainSubsetSCE[1:100, ],
                    reducedDimName = "PCA_counts"),
            "ggplot")
  expect_error(plotPCA(mouseBrainSubsetSCE[1:100, ],
                       reducedDimName = "missing"),
               "missing dimension not found. Run getPCA\\(\\) or set runPCA to TRUE.")
  expect_is(plotPCA(mouseBrainSubsetSCE[1:100, ],
                    reducedDimName = "missing", runPCA = TRUE),
            "ggplot")
  expect_error(plotPCA(mouseBrainSubsetSCE[1:100, ],
                       reducedDimName = "PCA_counts", pcX = "nope"),
               "pcX dimension nope is not in the reducedDim data")
  expect_error(plotPCA(mouseBrainSubsetSCE[1:100, ],
                       reducedDimName = "PCA_counts", pcY = "nope"),
               "pcY dimension nope is not in the reducedDim data")
})

test_that("getTSNE", {
  expect_error(getTSNE(mouseBrainSubsetSCE, "abcd"),
               "abcd not in the assay list")
  expect_error(getTSNE(mouseBrainSubsetSCE[1:100, ], "abcd"),
               "abcd not in the assay list")
  expect_is(getTSNE(mouseBrainSubsetSCE[1:100, ], "counts"),
            "SCtkExperiment")
  assay(mouseBrainSubsetSCE, "test") <- data.frame(assay(mouseBrainSubsetSCE,
                                                         "counts"))
  expect_error(getTSNE(mouseBrainSubsetSCE[1:100, ], "test"),
               "Input matrix test is not a matrix")
})

test_that("plotTSNE", {
  expect_is(plotTSNE(mouseBrainSubsetSCE[1:100, ],
                    reducedDimName = "TSNE_counts"),
            "ggplot")
  expect_error(plotTSNE(mouseBrainSubsetSCE[1:100, ],
                       reducedDimName = "missing"),
               "missing dimension not found. Run getTSNE\\(\\) or set runTSNE to TRUE.")
  expect_is(plotTSNE(mouseBrainSubsetSCE[1:100, ],
                    reducedDimName = "missing", runTSNE = TRUE),
            "ggplot")
})

test_that("getBiomarker", {
  expect_is(getBiomarker(mouseBrainSubsetSCE[1:100, ], "Tspan12"),
            "data.frame")
  expect_is(getBiomarker(mouseBrainSubsetSCE[1:100, ], "Tspan12",
                         binary = "Continuous"),
            "data.frame")
})

test_that("plotBiomarker", {
  expect_is(plotBiomarker(mouseBrainSubsetSCE[1:100, ],
                          c("Tspan12", "Tshz1"), binary = "Continuous",
                          reducedDimName="TSNE_counts"),
                "NULL")
  expect_error(plotBiomarker(mouseBrainSubsetSCE[1:100, ],
                          c("Tspan12"), binary = "Continuous",
                          reducedDimName="TSNE"),
            "Please supply a correct reducedDimName")
})
