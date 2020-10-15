#' @title Detecting contamination with DecontX.
#' @description A wrapper function for \link[celda]{decontX}. Identify
#'  potential contamination from experimental factors such as ambient RNA.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' Default NULL. \link[celda]{decontX} will be run on cells from each
#' sample separately.
#' @param useAssay  A string specifying which assay in the SCE to use. Default
#' 'counts'.
#' @param z Numeric or character vector. Cell cluster labels. If NULL,
#' PCA will be used to reduce the dimensionality of the dataset initially,
#' '\link[uwot]{umap}' from the 'uwot' package
#' will be used to further reduce the dataset to 2 dimenions and
#' the '\link[dbscan]{dbscan}' function from the 'dbscan' package
#' will be used to identify clusters of broad cell types. Default NULL.
#' @param maxIter Integer. Maximum iterations of the EM algorithm. Default 500.
#' @param convergence Numeric. The EM algorithm will be stopped if the maximum
#' difference in the contamination estimates between the previous and
#' current iterations is less than this. Default 0.001.
#' @param iterLogLik Integer. Calculate log likelihood every \code{iterLogLik}
#' iteration. Default 10.
#' @param delta Numeric Vector of length 2. Concentration parameters for
#' the Dirichlet prior for the contamination in each cell. The first element
#' is the prior for the native counts while the second element is the prior for
#' the contamination counts. These essentially act as pseudocounts for the
#' native and contamination in each cell. If \code{estimateDelta = TRUE},
#' this is only used to produce a random sample of proportions for an initial
#' value of contamination in each cell. Then
#' \code{\link[MCMCprecision]{fit_dirichlet}} is used to update
#' \code{delta} in each iteration.
#' If \code{estimateDelta = FALSE}, then \code{delta} is fixed with these
#' values for the entire inference procedure. Fixing \code{delta} and
#' setting a high number in the second element will force \code{decontX}
#' to be more aggressive and estimate higher levels of contamination at
#' the expense of potentially removing native expression.
#' Default \code{c(10, 10)}.
#' @param estimateDelta Boolean. Whether to update \code{delta} at each
#' iteration.
#' @param varGenes Integer. The number of variable genes to use in
#' dimensionality reduction before clustering. Variability is calcualted using
#' \code{\link[scran]{modelGeneVar}} function from the 'scran' package.
#' Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter
#' used in '\link[dbscan]{dbscan}' to estimate broad cell clusters.
#' Used only when z is not provided. Default 1.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  'decontX_Contamination' and 'decontX_Clusters' added to the
#'  \link[SummarizedExperiment]{colData} slot. Additionally, the
#' decontaminated counts will be added as an assay called 'decontXCounts'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDecontX(sce)
#' @export
runDecontX <- function(inSCE,
    sample = NULL,
    useAssay = "counts",
    z = NULL,
    maxIter = 500,
    delta = c(10, 10),
    estimateDelta = TRUE,
    convergence = 0.001,
    iterLogLik = 10,
    varGenes = 5000,
    dbscanEps = 1,
    seed = 12345,
    logfile = NULL,
    verbose = TRUE
) {
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  if(!is.null(sample)) {
    if(length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  }

  message(paste0(date(), " ... Running 'DecontX'"))

  inSCE <- celda::decontX(x = inSCE,
                          batch = sample,
                          assayName = useAssay,
                          z = z,
                          maxIter = maxIter,
                          delta = delta,
                          estimateDelta = estimateDelta,
                          convergence = convergence,
                          iterLogLik = iterLogLik,
                          varGenes = varGenes,
                          dbscanEps = dbscanEps,
                          seed = seed,
                          logfile = logfile,
                          verbose = verbose)

  #argsList <- argsList[!names(argsList) %in% ("...")]

  inSCE@metadata$runDecontX <- argsList[-1]
  inSCE@metadata$runDecontX$packageVersion <- utils::packageDescription("celda")$Version

  return(inSCE)
}

