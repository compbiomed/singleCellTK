context("check for conditions to perform enrichment using enrichR")

test_that("Should fail when the object supplied is not a sctk/ sctkexperiment object", {
  expect_error(enrichRSCE(matrix(1:100), "Cmtm5", NULL),
               "Please use a singleCellTK or a SCtkExperiment object")
})

test_that("Should fail when an incorrect gene name is provided", {
  expect_error(enrichRSCE(mouseBrainSubsetSCE, "mohammed", NULL),
               "Gene in gene list not found in input object.")
})

test_that("Should fail when an incorrect db name is provided", {
  expect_error(enrichRSCE(mouseBrainSubsetSCE, "Cmtm5", "somedb"),
  "database 'somedb' does not exist")
})

test_that("Should fail when the gene list is empty", {
  expect_error(enrichRSCE(mouseBrainSubsetSCE, NULL, NULL),
               "Please provide a gene list.")
})
