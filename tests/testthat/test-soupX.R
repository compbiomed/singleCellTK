library(singleCellTK)
context("Testing SoupX clustering")

set.seed(1)
nGenes <- 600
nCells <- 300
means <- rep(1, nGenes)
counts <- matrix(rpois(nGenes * nCells, 1), nrow = nGenes,
                 dimnames = list(paste0("g", seq_len(nGenes)),
                                 paste0("c", seq_len(nCells))))
groups <- rep(1:3, each = nCells / 3)
for (g in 1:3) {
  markers <- ((g - 1) * 50 + 1):(g * 50)
  counts[markers, groups == g] <- rpois(50 * sum(groups == g), 15)
}

test_that("quick clustering for SoupX recovers groups without deprecations", {
  expect_no_warning(cl <- .quickClusterRNA(counts),
                    class = "deprecatedWarning")
  expect_s3_class(cl, "factor")
  expect_length(cl, nCells)
  expect_true(all(table(cl) >= 100))
  ari <- bluster::pairwiseRand(groups, as.integer(cl), mode = "index")
  expect_gt(ari, 0.95)
})

test_that("quick clustering for SoupX stops on too few cells", {
  expect_error(.quickClusterRNA(counts[, 1:50]),
               "fewer cells than the minimum cluster size")
})

test_that("small clusters merge into the neighbor that best keeps modularity", {
  g <- igraph::make_graph(c(1, 2, 2, 3, 3, 1, 4, 5, 5, 6, 6, 4, 3, 7),
                          directed = FALSE)
  igraph::E(g)$weight <- 1
  merged <- .mergeSmallClusters(g, c(1, 1, 1, 2, 2, 2, 3), minSize = 2)
  expect_equal(merged, c(1L, 1L, 1L, 2L, 2L, 2L, 1L))
})
