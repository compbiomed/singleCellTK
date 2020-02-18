
.readDropEstFile <- function(dropEstdirs, counts_matrix_type,rdsFileName){
  dropEst_cell_counts <- paste0(dropEstdirs,'/',rdsFileName,'.rds')
  if (!file.exists(dropEst_cell_counts)){
    stop("DropEst output not found at location specified. Please check path provided and/or filename.")
  }
  dropEst_rds <- readRDS(dropEst_cell_counts) 
  
  return(dropEst_rds)
}

.constructColdata <- function(dropEst_rds,counts_matrix, counts_matrix_type){
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
  if (counts_matrix_type == 'filtered'){
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

.importDropEst <- function(dropEstdirs, counts_matrix_type, rdsFileName){
  ## Read DropEst RDS
  dropEst_rds <- .readDropEstFile(dropEstdirs,counts_matrix_type,rdsFileName)
  if (counts_matrix_type == 'filtered' && 'cm' %in% names(dropEst_rds)) {
    counts_matrix <- dropEst_rds$cm
  } else if (counts_matrix_type == 'raw' && 'cm_raw' %in% names(dropEst_rds)) {
    counts_matrix <- dropEst_rds$cm_raw
  } else {
    stop("No counts matrix found in the .rds provided! Exiting.")
  }
  ## Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix))
  
  ## Add SCE ColData. If using filtered counts matrix, colData is subset to include filtered cells.
  sce_coldata <- .constructColdata(dropEst_rds, counts_matrix, counts_matrix_type)
  if (dim(counts_matrix)[2] == dim(sce_coldata)[1]){
    SummarizedExperiment::colData(sce) <- sce_coldata
  } else {
    warning("Unable to add ColData to SCE. nCol of Counts Matrix not equal to nRow of ColData matrix.")
  }
  
  ## Add SCE metadata
  sce_metadata <- .extractMetadata(dropEst_rds)
  metadata(sce) <- sce_metadata
  
  return(sce)
}

#' @name importDropEst
#' @rdname importDropEst
#' @title Create a SingleCellExperiment Object from DropEst output 
#' @description Read the RDS file created by DropEst (https://github.com/hms-dbmi/dropEst) and
#' create a SingleCellExperiment object from either the raw or filtered counts matrix.
#' Additionally parse through the RDS to obtain appropriate feature annotations as 
#' SCE coldata, in addition to any metadata.
#' @param dropEstdirs Path to the folder containing the DropEst cell.counts.rds output
#' @param counts_matrix_type can be "filtered" or "raw". Default is "filtered"
#' @param rdsFileName File name prefix of the DropEst RDS output. default is "cell.counts"
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
#' # Example #1
#' Example DropEst outputs were downloaded from the DropEst Github 
#' (http://pklab.med.harvard.edu/viktor/dropest_paper/dropest_0.8.5.zip). 
#' To run the dropest import function with the example dataset, 
#' set the dropEstdirs variable to the example dropEst provided in SCTK as follows-
#' sce <- importDropEst(dropEstdirs = 'path/to/dropest/folder/', 
#'                      counts_matrix_type='filtered')
#' @export
importDropEst <- function(dropEstdirs = NULL, 
                          counts_matrix_type = 'filtered',
                          rdsFileName = 'cell.counts') {
  .importDropEst(dropEstdirs = dropEstdirs, 
                 counts_matrix_type = counts_matrix_type,
                 rdsFileName = rdsFileName)
}
  
  
### Test ####
#library("SingleCellExperiment")
#library("SummarizedExperiment")
#library("tidyverse")
#sce <- importDropEst(dropEstdirs = ../inst/extdata/DropEst_neurons_900_10x/', counts_matrix_type='filtered')
#print(names(metadata(sce)))


