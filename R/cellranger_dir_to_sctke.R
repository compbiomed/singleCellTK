#' Convert Cell Ranger outs Directory to a SingleCelltkExperiment
#'
#' @param path The path to the Cell Ranger outs directory. Defaults to the cwd.
#'
#' @return a SingleCelltkExperiment object to use in the Single Cell Toolkit.
#' raw data, PCA, tSNE, and clustering results are extracted and stored in the
#' object.
#' path <- "~/Dropbox/grad_school/johnson_lab/20170919_CellrangertoSCTK/JC-PBMC-HIGH/outs/"
#' @export
#'
cellranger_dir_to_sctke <- function(path=".") {
  if (!requireNamespace("cellrangerRkit", quietly = TRUE)) {
    stop("cellrangerRkit needed for this function to work. Please install it. ",
         "Directions are available here:\nhttps://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit",
         call. = FALSE)
  }
  
  if(length(list.files(path=file.path(path, "filtered_gene_bc_matrices"))) == 1){
    mat_dir <- list.files(path=file.path(path, "filtered_gene_bc_matrices"))
    barcode_file <- file.path(path, "filtered_gene_bc_matrices", mat_dir, "barcodes.tsv")
    genes_file <- file.path(path, "filtered_gene_bc_matrices", mat_dir, "genes.tsv")
    matrix_file <- file.path(path, "filtered_gene_bc_matrices", mat_dir, "matrix.mtx")
    if(file.exists(file.path(path, "metrics_summary.csv"))){
      metrics <- file.path(path, "metrics_summary.csv")
    } else {
      metrics <- NULL
    }
    incrdata <- cellrangerRkit::load_cellranger_matrix_from_files(mat_fn = matrix_file,
                                                                  gene_fn = genes_file,
                                                                  barcode_fn = barcode_file,
                                                                  summary_fn = metrics)
    inscedata <- createSCE(countfile = as.matrix(exprs(incrdata)),
                           annotfile = Biobase::pData(incrdata),
                           featurefile = Biobase::fData(incrdata),
                           inputdataframes = T)
    rm(incrdata)
  } else {
    stop("Error while parsing filtered_gene_bc_matrices, expected 1 sub directory.")
  }
  
  # if there is a pca folder
  if (dir.exists(file.path(path, "analysis", "pca"))){
    if(length(list.files(path=file.path(path, "analysis", "pca"))) > 1){
      warning("More than one pca components file found. Will consider first only.")
    }
    pca_dir <- list.files(pattern = "components", path=file.path(path, "analysis", "pca"))[1]
    pca_p <- read.csv(file.path(path, "analysis", "pca", pca_dir, "projection.csv"), row.names=1)
    pca_v <- read.csv(file.path(path, "analysis", "pca", pca_dir, "variance.csv"), row.names=1)
    
  }
  
  # if there is a tsne folder
  if (dir.exists(file.path(path, "analysis", "tsne"))){
    1+1
  }
  
  # if there is a clustering folder
  if (dir.exists(file.path(path, "analysis", "clustering"))){
    1+1
  }
  
  list.files(pattern = "components", path=file.path(path, "analysis", "pca"))
  return(inscedata)
}
