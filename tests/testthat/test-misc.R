context("misc functions")

test_that("Create SCTKE", {
  expect_is(createSCE(assayFile = assay(mouseBrainSubsetSCE[1:100, ]),
                      inputDataFrames = TRUE),
            "SCtkExperiment")
  expect_error(createSCE(assayFile = assay(mouseBrainSubsetSCE[1:100, ]),
                         annotFile = colData(mouseBrainSubsetSCE)[1:10, ],
                         inputDataFrames = TRUE),
               "Different number of samples in input matrix and annotations: annot: 10, counts: 30")
})

test_that("summarizeSCE", {
  ta <- summarizeSCE(mouseBrainSubsetSCE, sample = NULL)
  expect_is(ta, "data.frame")
})
