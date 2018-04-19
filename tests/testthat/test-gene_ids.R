context("gene id conversion tests")

test_that("gene id should fail if a bad database is given", {
  expect_error(convertGeneIDs(mouseBrainSubsetSCE, inSymbol = "SYMBOL",
                              outSymbol = "ENSEMBL", "this_should_fail"),
               "The database you want to use, this_should_fail, is not supported")
})

if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  library(org.Mm.eg.db)
  test_that("conversion of mouse brain subset should be successful", {
    expect_equal(ncol(convertGeneIDs(mouseBrainSubsetSCE,
                                    inSymbol = "SYMBOL", outSymbol = "SYMBOL",
                                    database = "org.Mm.eg.db")),
                 30)
    expect_equal(ncol(convertGeneIDs(mouseBrainSubsetSCE,
                                    inSymbol = "SYMBOL",
                                    outSymbol = "ENSEMBL",
                                    database =  "org.Mm.eg.db")),
                 30)
    expect_equal(ncol(convertGeneIDs(mouseBrainSubsetSCE,
                                    inSymbol = "SYMBOL",
                                    outSymbol = "ENTREZID", "org.Mm.eg.db")),
                 30)
  })
}
