#' @name importAlevin
#' @rdname importAlevin
#' @title Construct SCE object from Salmon-Alevin output
#' @param alevinDir Character. The output directory of salmon-Alevin pipeline.
#'  It should contain subfolder named 'alevin', which contains the count data which is stored
#'  in 'quants_mat.gz'. Default \code{NULL}.
#' @param sampleName Character. A user-defined sample name for the sample to be
#'  imported. The 'sampleName' will be appended to the begining of cell barcodes. Default is 'sample'.
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link{DelayedArray} object or not. Default \code{FALSE}.
#' @return A \code{SingleCellExperiment} object containing the count
#'  matrix, the feature annotations, and the cell annotation
#'  (which includes QC metrics stored in 'featureDump.txt').
#' @export

importAlevin <- function(
	alevinDir = NULL,
	sampleName = 'sample',
	delayedArray = FALSE) {

	matFile <- file.path(alevinDir, "alevin/quants_mat.gz")
	ma <- tximport::tximport(files = matFile, type = "alevin") ### require package fishpond
	if ('counts' %in% names(ma)) {
		stop("RNA count matrix not found in the alevin output!")
	}
	mat <- ma$counts
	if (delayedArray) {
		mat <- DelayedArray::DelayedArray(mat)
	}
	genes <- rownames(mat)
	cb <- .readBarcodes(file.path(alevinDir, 'alevin/featureDump.txt'),
						header = 'auto',
              			colname = "cell_barcode",
              			colClasses = "character")
	coln <- paste(sampleName, cb[[1]], sep = "_")

	sce <- SingleCellExperiment::SingleCellExperiment(
  											assays = list(counts = mat))
	SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(
											'feature_name' = genes,
                                   			row.names = genes)
	SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
											cb,
											column_name = coln,
											sample = sampleName,
											row.names = coln)

	return(sce)
}
