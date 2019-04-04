context("GSVA Testing")

test_that("GSVA Error Messages", {
  expect_error(gsvaPlot(plotType = "noplot"),
               "ERROR: Unsupported plot type")
  expect_error(gsvaPlot(gsvaData = data.frame(res = 1:100), plotType = "Violin"),
               "Too Many results for Violin Plot. Try Heatmap.")
  expect_error(gsvaPlot(gsvaData = data.frame(res = 1:40), plotType = "Violin"),
               "You must specify a condition for Violin plot")
})

utils::data(maits, package = "MAST")
utils::data(c2BroadSets, package = "GSVAdata")
maitslogtpm <- t(maits$expressionmat)
genesToSubset <- rownames(maitslogtpm)[which(rownames(maitslogtpm) %in%
                                       GSEABase::geneIds(c2BroadSets[["KEGG_PROTEASOME"]]))]
maitslogtpm <- maitslogtpm[rownames(maitslogtpm) %in% genesToSubset, ]
maitsfeatures <- maits$fdat[rownames(maits$fdat) %in% genesToSubset, ]
maitsSCE <- createSCE(assayFile = maitslogtpm,
                      annotFile = maits$cdat,
                      featureFile = maitsfeatures,
                      assayName = "logtpm",
                      inputDataFrames = TRUE,
                      createLogCounts = FALSE)
rowData(maitsSCE)$testbiomarker <- rep(1, nrow(maitsSCE))

test_that("GSVA Works", {
  expect_is(gsvaSCE(inSCE = maitsSCE, useAssay = "logtpm",
                    pathwaySource = "Manual Input",
                    pathwayNames = "testbiomarker", parallel.sz = 1),
            "matrix")
  expect_is(gsvaSCE(inSCE = maitsSCE, useAssay = "logtpm",
                    pathwaySource = "MSigDB c2 (Human, Entrez ID only)",
                    pathwayNames = "KEGG_PROTEASOME", parallel.sz = 1),
            "matrix")
  bigres <- gsvaSCE(inSCE = maitsSCE, useAssay = "logtpm",
                    pathwaySource = "MSigDB c2 (Human, Entrez ID only)",
                    pathwayNames = "ALL", parallel.sz = 1)
  expect_is(gsvaPlot(inSCE = maitsSCE, gsvaData = bigres[1:2, ],
                     plotType = "Violin", condition = "condition"),
            "ggplot")
  expect_is(gsvaPlot(inSCE = maitsSCE, gsvaData = bigres[1:2, ],
                     plotType = "Heatmap"),
            "HeatmapList")
})
