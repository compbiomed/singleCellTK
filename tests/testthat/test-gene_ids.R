context("gene id conversion tests")

test_that("gene id should fail if a bad database is given", {
  expect_error(convert_gene_ids(GSE60361_subset_sce, in_symbol = "SYMBOL",
                                out_symbol = "ENSEMBL", "this_should_fail"),
               "The database you want to use, this_should_fail, is not supported")
})

if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  library(org.Mm.eg.db)
  test_that("conversion of GSE60361 should be successful", {
    expect_equal(dim(convert_gene_ids(GSE60361_subset_sce, in_symbol = "SYMBOL",
                                      out_symbol = "SYMBOL",
                                      database = "org.Mm.eg.db")),
                 c(19972, 30))
    expect_equal(dim(convert_gene_ids(GSE60361_subset_sce, in_symbol = "SYMBOL",
                                      out_symbol = "ENSEMBL",
                                      database =  "org.Mm.eg.db")),
                 c(18153, 30))
    expect_equal(dim(convert_gene_ids(GSE60361_subset_sce, in_symbol = "SYMBOL",
                                      out_symbol = "ENTREZID", "org.Mm.eg.db")),
                 c(18735, 30))
  })
}
