
#' @title Find doublets/multiplets using \link[scds]{cxds}.
#' @description A wrapper function for \link[scds]{cxds}. Annotate
#'  doublets/multiplets using co-expression based approach. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[scds]{cxds} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{cxds}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{cxds_score} and optionally \emph{cxds_call}.
#'  Please refer to the documentation of \link[scds]{cxds} for details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runCxds(sce_chcl)
#' @export
#' @import scds
runCxds <- function(inSCE,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'cxds'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            cxds_score = numeric(ncol(inSCE)),
            cxds_call = logical(ncol(inSCE)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            cxds_score = numeric(ncol(inSCE)))
    }

    ## Loop through each sample and run cxds
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- inSCE[, sceSampleInd]

        counts(sceSample) <- .convertToMatrix(counts(sceSample))
        
        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed, scds::cxds(sce = sceSample, ntop = nGene, ...))}, silent = TRUE)
          nGene <- nGene - 100
        }  
        
        if (!inherits(result, "try-error") & !is.null(result)) {
          if ("cxds_call" %in% colnames(SummarizedExperiment::colData(result))) {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("cxds_score", "cxds_call")]
          } else {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("cxds_score")]
          }
        } else {
          output[sceSampleInd, ] <- NA
          warning(paste0("'cxds' from package 'scds' did not complete successfully for sample", samples[i]))
        }            
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(inSCE) = cbind(colData(inSCE), output)

    return(inSCE)
}


#' @title Find doublets/multiplets using \link[scds]{bcds}.
#' @description A wrapper function for \link[scds]{bcds}. Annotate
#'  doublets/multiplets using a binary classification approach to discriminate
#'  artificial doublets from original data. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[scds]{bcds} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{bcds}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{bcds} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{bcds_score} and optionally \emph{bcds_call}.
#'  Please refer to the documentation of \link[scds]{bcds} for details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runBcds(sce_chcl)
#' @export
#' @import scds
runBcds <- function(inSCE,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'bcds'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            bcds_score = numeric(ncol(inSCE)),
            bcds_call = logical(ncol(inSCE)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            bcds_score = numeric(ncol(inSCE)))
    }

    ## Loop through each sample and run bcds
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- inSCE[, sceSampleInd]

        counts(sceSample) <- .convertToMatrix(counts(sceSample))
        
        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed, scds::bcds(sce = sceSample, ntop = nGene, ...))}, silent = TRUE)
          nGene <- nGene - 100
        }  

        if (!inherits(result, "try-error") & !is.null(result)) {
          if ("bcds_call" %in% colnames(SummarizedExperiment::colData(result))) {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("bcds_score", "bcds_call")]
          } else {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("bcds_score")]
          }
        } else {
          output[sceSampleInd, ] <- NA
          warning(paste0("'bcds' from package 'scds' did not complete successfully for sample", samples[i]))
        }  
 
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(inSCE) = cbind(colData(inSCE), output)

    return(inSCE)
}


#' @title Find doublets/multiplets using \link[scds]{cxds_bcds_hybrid}.
#' @description A wrapper function for \link[scds]{cxds_bcds_hybrid}. Annotate
#'  doublets/multiplets using a binary classification approach to discriminate
#'  artificial doublets from original data. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[scds]{cxds_bcds_hybrid} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{cxds_bcds_hybrid}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds_bcds_hybrid} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{hybrid_score} and optionally \emph{hybrid_call}.
#'  Please refer to the documentation of \link[scds]{cxds_bcds_hybrid} for
#'  details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runCxdsBcdsHybrid(sce_chcl)
#' @export
#' @import scds
runCxdsBcdsHybrid <- function(inSCE,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'cxds_bcds_hybrid'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            hybrid_score = numeric(ncol(inSCE)),
            hybrid_call = logical(ncol(inSCE)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(inSCE),
            hybrid_score = numeric(ncol(inSCE)))
    }

    ## Loop through each sample and run cxds_bcds_hybrid
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- inSCE[, sceSampleInd]

        counts(sceSample) <- .convertToMatrix(counts(sceSample))

        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed, scds::cxds_bcds_hybrid(sce = sceSample, cxdsArgs=list(ntop = nGene), bcdsArgs=list(ntop = nGene)))}, silent = TRUE)
          nGene <- nGene - 100
        }  

        if (!inherits(result, "try-error") & !is.null(result)) {
          if ("hybrid_call" %in% colnames(SummarizedExperiment::colData(result))) {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("hybrid_score", "hybrid_call")]
          } else {
              output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                  c("hybrid_score")]
          }
        } else {
          output[sceSampleInd, ] <- NA
          warning(paste0("'cxds_bcds_hybrid' from package 'scds' did not complete successfully for sample", samples[i]))
        }   
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(inSCE) = cbind(colData(inSCE), output)

    return(inSCE)
}
