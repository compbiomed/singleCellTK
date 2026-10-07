library(singleCellTK)
context("Testing per-cell QC metrics")

set.seed(1)
counts <- matrix(rpois(40 * 6, 2), nrow = 40,
                 dimnames = list(paste0("g", 1:40), paste0("c", 1:6)))
counts[, 2] <- 0
sce <- SingleCellExperiment(list(counts = counts))
spikes <- matrix(rpois(5 * 6, 3), nrow = 5,
                 dimnames = list(paste0("s", 1:5), colnames(counts)))
altExp(sce, "ERCC") <- SingleCellExperiment(list(counts = spikes))
mito <- c("g1", "g2", "g3")

test_that("runPerCellQC metrics match base R without deprecations", {
  expect_no_warning(
    res <- runPerCellQC(sce, geneSetList = list(mito = mito),
                        geneSetListLocation = "rownames",
                        percent_top = c(5, 50), use_altexps = TRUE),
    class = "deprecatedWarning"
  )
  cd <- colData(res)
  libSize <- colSums(counts)
  top5 <- apply(counts, 2, function(v) sum(sort(v, decreasing = TRUE)[1:5]))
  total <- libSize + colSums(spikes)
  expect_equal(cd$sum, unname(libSize))
  expect_equal(cd$detected, unname(colSums(counts > 0)))
  expect_equal(cd$percent.top_5, unname(top5 / libSize * 100))
  expect_equal(cd$percent.top_50, unname(ifelse(libSize > 0, 100, NaN)))
  expect_equal(cd$mito_sum, unname(colSums(counts[mito, ])))
  expect_equal(cd$mito_detected, unname(colSums(counts[mito, ] > 0)))
  expect_equal(cd$mito_percent, unname(colSums(counts[mito, ]) / libSize * 100))
  expect_equal(cd$altexps_ERCC_sum, unname(colSums(spikes)))
  expect_equal(cd$altexps_ERCC_percent, unname(colSums(spikes) / total * 100))
  expect_equal(cd$total, unname(total))
  expectedNames <- c("sum", "detected", "percent.top_5", "percent.top_50",
                     "mito_sum", "mito_detected", "mito_percent",
                     "altexps_ERCC_sum", "altexps_ERCC_detected",
                     "altexps_ERCC_percent", "total")
  expect_equal(setdiff(names(cd), names(colData(sce))), expectedNames)
})

test_that("runPerCellQC applies detectionLimit and flatten = FALSE", {
  res <- suppressWarnings(runPerCellQC(sce, geneSetList = list(mito = mito),
                                       geneSetListLocation = "rownames",
                                       percent_top = 5, detectionLimit = 2,
                                       flatten = FALSE))
  cd <- colData(res)
  expect_equal(cd$detected, unname(colSums(counts > 2)))
  expect_equal(cd$subsets$mito$detected,
               unname(colSums(counts[mito, ] > 2)))
  expect_equal(colnames(cd$percent.top), "5")
})

test_that("sampleSummaryStats adds sum and detected without deprecations", {
  expect_no_warning(res <- sampleSummaryStats(sce),
                    class = "deprecatedWarning")
  expect_equal(res$sum, unname(colSums(counts)))
  expect_equal(res$detected, unname(colSums(counts > 0)))
})

test_that("sampleSummaryStats uses only altExps that have counts", {
  sceAlt <- sce
  scaled <- matrix(rnorm(30), nrow = 5,
                   dimnames = list(NULL, colnames(counts)))
  altExp(sceAlt, "scaled") <- SingleCellExperiment(list(scaledata = scaled))
  res <- sampleSummaryStats(sceAlt)
  expect_true("altexps_ERCC_sum" %in% names(colData(res)))
  expect_false(any(grepl("altexps_scaled", names(colData(res)))))
  expect_equal(res$total, unname(colSums(counts) + colSums(spikes)))
})

test_that("QC metrics match base R for chunked, delayed, and decimal input", {
  decimal <- counts * 0.5
  decimal[3, 4] <- 0.25
  delayed <- DelayedArray::DelayedArray(decimal)
  sceDelayed <- SingleCellExperiment(list(counts = delayed))
  qc <- .perCellQCMetrics(sceDelayed, subsets = list(mito = mito),
                          percentTop = c(5, 50), chunkSize = 2,
                          numThreads = 2)
  libSize <- colSums(decimal)
  top5 <- apply(decimal, 2, function(v) sum(sort(v, decreasing = TRUE)[1:5]))
  expect_equal(qc$sum, unname(libSize))
  expect_equal(qc$detected, unname(colSums(decimal > 0)))
  expect_equal(qc$percent.top_5, unname(top5 / libSize * 100))
  expect_equal(qc$subsets_mito_sum, unname(colSums(decimal[mito, ])))
})
