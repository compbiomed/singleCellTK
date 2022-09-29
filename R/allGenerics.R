#' @title Stores and returns table of SCTK QC outputs to metadata.
#' @rdname getSampleSummaryStatsTable
#' @description  Stores and returns table of QC metrics generated from
#'  QC algorithms within the metadata slot of the SingleCellExperiment object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' \link{assay} data and/or \link{colData} data. Required.
#' @param statsName A \code{character} value indicating the slot
#' that stores the stats table within the metadata of the
#' SingleCellExperiment object. Required.
#' @param ... Other arguments passed to the function. 
#' @return For \code{getSampleSummaryStatsTable}, A matrix/array object. 
#' Contains a summary table for QC statistics generated from SingleCellTK. For
#' \code{setSampleSummaryStatsTable<-}, A SingleCellExperiment object where the 
#' summary table is updated in the \code{metadata} slot.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- sampleSummaryStats(sce, simple = TRUE, statsName = "qc_table")
#' getSampleSummaryStatsTable(sce, statsName = "qc_table")
#' @export
setGeneric("getSampleSummaryStatsTable", function(inSCE, statsName, ...) standardGeneric("getSampleSummaryStatsTable"))

#' @rdname getSampleSummaryStatsTable
setGeneric("setSampleSummaryStatsTable<-", function(inSCE, statsName, ..., value) standardGeneric("setSampleSummaryStatsTable<-"))

#' @title Lists the table of SCTK QC outputs stored within the metadata.
#' @rdname listSampleSummaryStatsTables
#' @description  Returns a character vector of the tables
#' within the metadata slot of the SingleCellExperiment object.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object with saved
#' table within the \link{metadata} data. Required.
#' @param ... Other arguments passed to the function. 
#' @return A character vector. Contains a list of summary tables
#' within the SingleCellExperiment object.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- sampleSummaryStats(sce, simple = TRUE, statsName = "qc_table")
#' listSampleSummaryStatsTables(sce)
#' @export
setGeneric("listSampleSummaryStatsTables", function(inSCE, ...) standardGeneric("listSampleSummaryStatsTables"))