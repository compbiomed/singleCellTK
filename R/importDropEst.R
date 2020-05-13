
.readDropEstFile <- function(sampleDir, dataType,rdsFileName){
  dropEst_cell_counts <- file.path(sampleDir, paste(rdsFileName, '.rds', sep=''))
  if (!file.exists(dropEst_cell_counts)){
    stop("DropEst output not found at location specified. Please check path provided and/or filename.")
  }
  dropEst_rds <- readRDS(dropEst_cell_counts)

  return(dropEst_rds)
}

.constructColdata <- function(dropEst_rds,counts_matrix, dataType){
  coldata_fields <- c("mean_reads_per_umi","aligned_reads_per_cell","aligned_umis_per_cell","requested_umis_per_cb","requested_reads_per_cb")
  coldata_df <-  list()
  for (field in coldata_fields){
    if (field %in% names(dropEst_rds)){
      coldata_field_df <- data.frame(as.matrix(dropEst_rds[[field]]))
      names(coldata_field_df)[1] <- field
      coldata_field_df$cell <- row.names(coldata_field_df)

      coldata_df[[field]] <- coldata_field_df
    }}
  coldata_df_merged <- Reduce(function(x, y) merge(x, y, all=TRUE,by="cell"), coldata_df)
  row.names(coldata_df_merged) <- coldata_df_merged$cell
  coldata_df_merged <- S4Vectors::DataFrame(as.matrix(coldata_df_merged))
  if (dataType == 'filtered'){
    coldata_df_merged <- coldata_df_merged[coldata_df_merged$cell %in% colnames(counts_matrix),]
  }
  return(coldata_df_merged)
}

.extractMetadata <- function(dropEst_rds){
  metadata_fields <- c("saturation_info","merge_targets","reads_per_umi_per_cell")
  metadata <- c()
  for (md in metadata_fields){
    if (md %in% names(dropEst_rds)){
      metadata[[md]] <- dropEst_rds[[md]]
    }}
  return(metadata)
}

.importDropEstSample <- function(sampleDir = './',
                                 dataType,
                                 rdsFileName,
                                 sampleName = 'sample',
                                 delayedArray = FALSE){
  ## Read DropEst RDS
  dropEst_rds <- .readDropEstFile(sampleDir,dataType,rdsFileName)
  if (dataType == 'filtered' && 'cm' %in% names(dropEst_rds)) {
    counts_matrix <- dropEst_rds$cm
  } else if (dataType == 'raw' && 'cm_raw' %in% names(dropEst_rds)) {
    counts_matrix <- dropEst_rds$cm_raw
  } else {
    stop("No counts matrix found in the .rds provided! Exiting.")
  }

  if (isTRUE(delayedArray)) {
    counts_matrix <- DelayedArray::DelayedArray(counts_matrix)
    }
  ## Create SingleCellExperiment object
  ## Add SCE ColData. If using filtered counts matrix, colData is subset to include filtered cells.
  ## append sample name to cells in SCE
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix))
  colnames(sce) <- paste0(sampleName,"_",colnames(sce))
  sce_coldata <- .constructColdata(dropEst_rds, counts_matrix, dataType)
  row.names(sce_coldata) <- paste0(sampleName,"_",row.names(sce_coldata))

  if (dim(counts_matrix)[2] == dim(sce_coldata)[1]){
    SummarizedExperiment::colData(sce) <- sce_coldata
  } else {
    warning("Unable to add ColData to SCE. nCol of Counts Matrix not equal to nRow of ColData matrix.")
  }

  ## Add SCE metadata
  sce_metadata <- .extractMetadata(dropEst_rds)
  sce@metadata$dropEst <- sce_metadata
  
  return(sce)
}

#' @name importDropEst
#' @rdname importDropEst
#' @title Create a SingleCellExperiment Object from DropEst output
#' @description imports the RDS file created by DropEst (https://github.com/hms-dbmi/dropEst) and
#' create a SingleCellExperiment object from either the raw or filtered counts matrix.
#' Additionally parse through the RDS to obtain appropriate feature annotations as
#' SCE coldata, in addition to any metadata.
#' @param sampleDirs  A path to the directory containing the data files. Default "./".
#' @param sampleNames A User-defined sample name. This will be prepended to all cell barcode IDs.
#'  Default "sample".
#' @param dataType can be "filtered" or "raw". Default \code{"filtered"}.
#' @param rdsFileName File name prefix of the DropEst RDS output. default is "cell.counts"
#' @param delayedArray Boolean. Whether to read the expression matrix as
#'  \link[DelayedArray]{DelayedArray} object or not. Default \code{TRUE}.
#' @details
#' \code{importDropEst} expects either raw counts matrix stored as "cm_raw" or filtered
#' counts matrix stored as "cm" in the DropEst rds output.
#' ColData is obtained from the DropEst corresponding to "mean_reads_per_umi","aligned_reads_per_cell",
#' "aligned_umis_per_cell","requested_umis_per_cb","requested_reads_per_cb"
#' If using filtered counts matrix, the colData dataframe is
#' subset to contain features from the filtered counts matrix alone.
#' If any annotations of ("saturation_info","merge_targets","reads_per_umi_per_cell") are
#' found in the DropEst rds, they will be added to the SCE metadata field
#' @return A \code{SingleCellExperiment} object containing the count matrix,
#'  the feature annotations from DropEst as ColData, and any metadata from DropEst
#' @examples
#' # Example results were generated as per instructions from the developers of dropEst described in
#' # https://github.com/hms-dbmi/dropEst/blob/master/examples/EXAMPLES.md
#' sce <- importDropEst(sampleDirs = system.file("extdata/dropEst_scg71", package = "singleCellTK"),
#'                      sampleNames = 'scg71')

#' @export
importDropEst <- function(sampleDirs = NULL,
                          dataType = c('filtered','raw'),
                          rdsFileName = 'cell.counts',
                          sampleNames = NULL,
                          delayedArray = TRUE) {
  dataType <- match.arg(dataType)

  if (length(sampleDirs)!=length(sampleNames)){
    stop("Please provide sample names for all input directories")
  }

  res <- vector("list", length = length(sampleDirs))

  for (i in seq_along(sampleDirs)){
    scei <- .importDropEstSample(sampleDir = sampleDirs[[i]],
                         sampleName = sampleNames[[i]],
                         dataType = dataType,
                         rdsFileName = rdsFileName,
                         delayedArray = delayedArray)
    res[[i]] <- scei
  }
  sce <- do.call(SingleCellExperiment::cbind, res)
  return(sce)
}





