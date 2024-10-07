#' @title Find doublets using \code{scrublet}.
#' @description A wrapper function that calls \code{scrub_doublets} from python
#'  module \code{scrublet}. Simulates doublets from the observed data and uses
#'  a k-nearest-neighbor classifier to calculate a continuous
#'  \code{scrublet_score} (between 0 and 1) for each transcriptome. The score
#'  is automatically thresholded to generate \code{scrublet_call}, a boolean
#'  array that is \code{TRUE} for predicted doublets and \code{FALSE}
#'  otherwise.
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param sample Character vector or colData variable name. Indicates which 
#' sample each cell belongs to. Default \code{NULL}.
#' @param useAssay  A string specifying which assay in the SCE to use. Default 
#' \code{"counts"}.
#' @param simDoubletRatio Numeric. Number of doublets to simulate relative to
#' the number of observed transcriptomes. Default \code{2.0}.
#' @param nNeighbors Integer. Number of neighbors used to construct the KNN
#' graph of observed transcriptomes and simulated doublets. If \code{NULL},
#' this is set to \code{round(0.5 * sqrt(n_cells))}. Default \code{NULL}.
#' @param minDist Float Determines how tightly UMAP packs points together. If 
#' \code{NULL}, this is set to \code{0.1}. Default \code{NULL}.
#' @param expectedDoubletRate The estimated doublet rate for the experiment.
#' Default \code{0.1}.
#' @param stdevDoubletRate Uncertainty in the expected doublet rate. Default 
#' \code{0.02}.
#' @param syntheticDoubletUmiSubsampling Numeric. Rate for sampling UMIs when
#' creating synthetic doublets. If \code{1.0}, each doublet is created by simply
#' adding the UMIs from two randomly sampled observed transcriptomes. For
#' values less than 1, the UMI counts are added and then randomly sampled at
#' the specified rate. Defuault \code{1.0}.
#' @param useApproxNeighbors Boolean. Use approximate nearest neighbor method
#' (annoy) for the KNN classifier. Default \code{TRUE}.
#' @param distanceMetric Character. Distance metric used when finding nearest
#' neighbors. See detail. Default \code{"euclidean"}. 
#' @param getDoubletNeighborParents Boolean. If \code{TRUE}, return the
#' parent transcriptomes that generated the doublet neighbors of each
#' observed transcriptome. This information can be used to infer the cell states
#' that generated a given doublet state. Default \code{FALSE}.
#' @param minCounts Numeric. Used for gene filtering prior to PCA. Genes
#' expressed at fewer than \code{minCounts} in fewer than \code{minCells} are
#' excluded. Default \code{3}.
#' @param minCells Integer. Used for gene filtering prior to PCA. Genes
#' expressed at fewer than \code{minCounts} in fewer than \code{minCells} are 
#' excluded. Default \code{3}.
#' @param minGeneVariabilityPctl Numeric. Used for gene filtering prior to
#' PCA. Keep the most highly variable genes (in the top
#' \code{minGeneVariabilityPctl} percentile), as measured by the v-statistic
#' (Klein et al., Cell 2015). Default \code{85}.
#' @param logTransform Boolean. If \code{TRUE}, log-transform the counts matrix
#' (log1p(TPM)). \code{sklearn.decomposition.TruncatedSVD} will be used for
#' dimensionality reduction, unless \code{meanCenter} is \code{TRUE}. Default 
#' \code{FALSE}.
#' @param meanCenter If \code{TRUE}, center the data such that each gene has a
#' mean of \code{0}. \code{sklearn.decomposition.PCA} will be used for
#' dimensionality reduction. Default \code{TRUE}.
#' @param normalizeVariance Boolean. If \code{TRUE}, normalize the data such
#' that each gene has a variance of 1.
#' \code{sklearn.decomposition.TruncatedSVD} will be used for dimensionality
#' reduction, unless \code{meanCenter} is \code{TRUE}. Default \code{TRUE}.
#' @param nPrinComps Integer. Number of principal components used to embed
#' the transcriptomes prior to k-nearest-neighbor graph construction.
#' Default \code{30}.
#' @param tsneAngle Float. Determines angular size of a distant node as measured
#' from a point in the t-SNE plot. If \code{NULL}, it is set to \code{0.5}. 
#' Default \code{NULL}.
#' @param tsnePerplexity Integer. The number of nearest neighbors that is used
#' in other manifold learning algorithms. If \code{NULL}, it is set to 30. 
#' Default \code{NULL}.
#' @param verbose Boolean. If \code{TRUE}, print progress updates. Default
#' \code{TRUE}.
#' @param seed Seed for the random number generator, can be set to \code{NULL}. 
#' Default \code{12345}.
#' @details For the list of valid values for \code{distanceMetric}, see the 
#' documentation for 
#' \href{https://github.com/spotify/annoy#full-python-api}{annoy} (if 
#' \code{useApproxNeighbors} is \code{TRUE}) or 
#' \href{https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.NearestNeighbors.html}{sklearn.neighbors.NearestNeighbors}
#' (if \code{useApproxNeighbors} is  \code{FALSE}). 
#' @return A \linkS4class{SingleCellExperiment} object with
#' \code{scrub_doublets} output appended to the \link{colData} slot. The columns
#' include \code{scrublet_score} and \code{scrublet_call}.
#' @seealso \code{\link{plotScrubletResults}}, \code{\link{runCellQC}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runScrublet(sce)
#' }
#' @export
#' @importFrom reticulate py_module_available py_set_seed import
#' @importFrom SummarizedExperiment colData colData<- assayNames assayNames<-
#' @importFrom SingleCellExperiment reducedDim<- reducedDims
#' @importFrom S4Vectors metadata metadata<-
runScrublet <- function(inSCE,
                        sample = NULL,
                        useAssay = "counts",
                        simDoubletRatio = 2.0,
                        nNeighbors = NULL,
                        minDist = NULL,
                        expectedDoubletRate = 0.1,
                        stdevDoubletRate = 0.02,
                        syntheticDoubletUmiSubsampling = 1.0,
                        useApproxNeighbors = TRUE,
                        distanceMetric = "euclidean",
                        getDoubletNeighborParents = FALSE,
                        minCounts = 3,
                        minCells = 3L,
                        minGeneVariabilityPctl = 85,
                        logTransform = FALSE,
                        meanCenter = TRUE,
                        normalizeVariance = TRUE,
                        nPrinComps = 30L,
                        tsneAngle = NULL,
                        tsnePerplexity = NULL,
                        verbose = TRUE,
                        seed = 12345) {
  
  if (!reticulate::py_module_available(module = "scrublet")) {
    warning("Cannot find python module 'scrublet', please install Conda and",
            " run sctkPythonInstallConda() or run sctkPythonInstallVirtualEnv().",
            "If one of these have been previously run to install the modules,",
            "make sure to run selectSCTKConda() or selectSCTKVirtualEnvironment(),",
            " respectively, if R has been restarted since the module installation.",
            " Alternatively, Scrublet can be installed on the local machine",
            "with pip (e.g. pip install scrublet) and then the 'use_python()'",
            " function from the 'reticulate' package can be used to select the",
            " correct Python environment.")
    return(inSCE)
  }
  
  if (!is.null(seed)) {
    reticulate::py_set_seed(seed = seed)
  }
  
  tempSCE <- inSCE
  #assayNames(inSCE)[which(useAssay %in% assayNames(inSCE))] <- "counts"
  #useAssay <- "counts"
  
  p <- paste0(date(), " ... Running 'scrublet'")
  message(p)
  
  ##  Getting current arguments values
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  argsList <- argsList[!names(argsList) %in% c("inSCE")]
  argsList$packageVersion <- pkg_resources$"require"("scrublet")[[1]]
  
  sample <- .manageCellVar(inSCE, var = sample)
  if (is.null(sample)) {
    sample <- rep(1, ncol(inSCE))
  }
  ## Define result matrix for all samples
  output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
                                 scrublet_score = numeric(ncol(inSCE)),
                                 scrublet_call = logical(ncol(inSCE)))
  
  ## Loop through each sample and run scrublet

  # try/catch block for when Scrublet fails (which can be often)

  # boolean flag variable for if the try/catch works
  successful = FALSE

  result <- tryCatch(
    {
        samples <- unique(sample)
        umapDims <- matrix(ncol = 2,
                         nrow = ncol(inSCE))
        rownames(umapDims) = colnames(inSCE)
        colnames(umapDims) = c("UMAP_1", "UMAP_2")
        
        tsneDims <- matrix(ncol = 2,
                          nrow = ncol(inSCE))
        rownames(tsneDims) = colnames(inSCE)
        colnames(tsneDims) = c("TSNE_1", "TSNE_2")
        
        for (s in samples) {
          sceSampleInd <- sample == s
          sceSample <- inSCE[, sceSampleInd]
          
          mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
          mat <- .convertToMatrix(mat)
          
          scr <- scrublet$Scrublet(counts_matrix = t(mat),
                                  sim_doublet_ratio = simDoubletRatio,
                                  n_neighbors = nNeighbors,
                                  expected_doublet_rate = expectedDoubletRate,
                                  stdev_doublet_rate = stdevDoubletRate)
          
          result <- scr$scrub_doublets(
            synthetic_doublet_umi_subsampling = syntheticDoubletUmiSubsampling,
            use_approx_neighbors = useApproxNeighbors,
            distance_metric = distanceMetric,
            get_doublet_neighbor_parents = getDoubletNeighborParents,
            min_counts = minCounts,
            min_cells = as.integer(minCells),
            min_gene_variability_pctl = minGeneVariabilityPctl,
            log_transform = logTransform,
            mean_center = meanCenter,
            normalize_variance = normalizeVariance,
            n_prin_comps = as.integer(nPrinComps),
            verbose = verbose)
          
          output[sceSampleInd, "scrublet_score"] <- result[[1]]
          output[sceSampleInd, "scrublet_call"] <- result[[2]]
          
          ## Extract UMAP and TSNE coordinates
          if (is.null(nNeighbors) && is.null(minDist)) {
            umap_coordinates <- scrublet$get_umap(scr$manifold_obs_)
          } else if (!is.null(nNeighbors) && !is.null(minDist)) {
            umap_coordinates <- scrublet$get_umap(
              scr$manifold_obs_,
              n_neighbors = as.integer(nNeighbors),
              min_dist = minDist
            )
          } else {
            warning("`nNeighbors` and `minDist` has to be specified or set to NULL",
                    " at the same time. Setting both to NULL.")
            umap_coordinates <- scrublet$get_umap(scr$manifold_obs_)
          }
          umapDims[sceSampleInd, ] <- umap_coordinates
          
          if (is.null(tsneAngle) && is.null(tsnePerplexity)) {
            tsne_coordinates <- scrublet$get_tsne(scr$manifold_obs_)
          } else if (!is.null(tsneAngle) && !is.null(tsnePerplexity)) {
            tsne_coordinates <- scrublet$get_tsne(
              scr$manifold_obs_,
              angle = tsneAngle,
              perplexity = as.integer(tsnePerplexity)
            )
          } else {
            warning("`tsneAngle` and `tsnePerplexity` has to be specified or set ", 
                    "to NULL at the same time. Setting both to NULL.")
            tsne_coordinates <- scrublet$get_tsne(scr$manifold_obs_)
          }
          tsneDims[sceSampleInd, ] <- tsne_coordinates
          if (!identical(samples, 1)) {
            metadata(inSCE)$sctk$runScrublet[[s]] <- argsList
          }
        }
        if (identical(samples, 1)) {
          metadata(inSCE)$sctk$runScrublet$all_cells <- argsList
        }
        colData(inSCE)$scrublet_score <- NULL
        colData(inSCE)$scrublet_call <- NULL
        
        colData(inSCE) = cbind(colData(inSCE), output)
        
        ## convert doublet call from TRUE/FALSE to Singlet/Doublet
        inSCE$scrublet_call <- as.factor(inSCE$scrublet_call)
        
        levels(inSCE$scrublet_call) <- list(Singlet = "FALSE", Doublet = "TRUE")
        
        reducedDim(inSCE,'scrublet_TSNE') <- tsneDims
        reducedDim(inSCE,'scrublet_UMAP') <- umapDims
        
        colData(tempSCE) <- colData(inSCE)
        metadata(tempSCE) <- metadata(inSCE)
        reducedDims(tempSCE) <- reducedDims(inSCE)
        
        successful = TRUE
        return(tempSCE)   
    },
    error=function(cond) {
        p <- paste0(date(), " ... Scrublet did not complete successfully; Returning SCE without changes. Scrublet error:")
        message(p)
        message(cond)
        return(inSCE)
    },
    finally={
        if (isTRUE(successful)) {
          return(tempSCE)
        }
        else {
          return(inSCE)
        }
    }
    # if (inherits(error, "try-error")) {
  #   warning("Scrublet did not complete successfully. Returning SCE without",
  #           " making any changes. Error given by Scrublet: \n\n", error)
  #   return(inSCE)
  # }
  )
  return(result)
}
