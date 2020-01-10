
#' @title Find doublets/multiplets using \link[scds]{cxds}.
#' @description A wrapper function for \link[scds]{cxds}. Annotate
#'  doublets/multiplets using co-expression based approach. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{cxds}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds} output appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
#'  \emph{cxds_score} and optionally \emph{cxds_call}.
#'  Please refer to the documentation of \link[scds]{cxds} for details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runCxds(sce_chcl)
#' @export
#' @import scds
runCxds <- function(sce,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(sce)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'sce'")
        }
    } else {
        sample <- rep(1, ncol(sce))
    }

    message(paste0(date(), " ... Running 'cxds'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            cxds_score = numeric(ncol(sce)),
            cxds_call = logical(ncol(sce)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            cxds_score = numeric(ncol(sce)))
    }

    ## Loop through each sample and run cxds
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- sce[, sceSampleInd]

        counts(sceSample) <- as(counts(sceSample), "dgCMatrix")
        result <- withr::with_seed(seed, scds::cxds(sce = sceSample, ...))

        if ("cxds_call" %in% colnames(SummarizedExperiment::colData(result))) {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("cxds_score", "cxds_call")]
        } else {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("cxds_score")]
        }
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(sce) = cbind(colData(sce), output)

    return(sce)
}


#' @title Find doublets/multiplets using \link[scds]{bcds}.
#' @description A wrapper function for \link[scds]{bcds}. Annotate
#'  doublets/multiplets using a binary classification approach to discriminate
#'  artificial doublets from original data. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{bcds}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{bcds} output appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
#'  \emph{bcds_score} and optionally \emph{bcds_call}.
#'  Please refer to the documentation of \link[scds]{bcds} for details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runBcds(sce_chcl)
#' @export
#' @import scds
runBcds <- function(sce,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(sce)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'sce'")
        }
    } else {
        sample <- rep(1, ncol(sce))
    }

    message(paste0(date(), " ... Running 'bcds'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            bcds_score = numeric(ncol(sce)),
            bcds_call = logical(ncol(sce)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            bcds_score = numeric(ncol(sce)))
    }

    ## Loop through each sample and run bcds
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- sce[, sceSampleInd]

        counts(sceSample) <- as(counts(sceSample), "dgCMatrix")
        result <- withr::with_seed(seed, scds::bcds(sce = sceSample, ...))

        if ("bcds_call" %in% colnames(SummarizedExperiment::colData(result))) {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("bcds_score", "bcds_call")]
        } else {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("bcds_score")]
        }
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(sce) = cbind(colData(sce), output)

    return(sce)
}


#' @title Find doublets/multiplets using \link[scds]{cxds_bcds_hybrid}.
#' @description A wrapper function for \link[scds]{cxds_bcds_hybrid}. Annotate
#'  doublets/multiplets using a binary classification approach to discriminate
#'  artificial doublets from original data. Generate a doublet
#'  score for each cell. Infer doublets if \code{estNdbl} is \code{TRUE}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#'  \link[DropletUtils]{emptyDrops} will be run on cells from each sample
#'  separately. If NULL, then all cells will be processed together.
#'  Default NULL.
#' @param seed Seed for the random number generator. Default 12345.
#' @param ... Additional arguments passed to \link[scds]{cxds_bcds_hybrid}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds_bcds_hybrid} output appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
#'  \emph{hybrid_score} and optionally \emph{hybrid_call}.
#'  Please refer to the documentation of \link[scds]{cxds_bcds_hybrid} for
#'  details.
#' @examples
#' data(sce_chcl, package = "scds")
#' sce <- runCxdsBcdsHybrid(sce_chcl)
#' @export
#' @import scds
runCxdsBcdsHybrid <- function(sce,
    sample = NULL,
    seed = 12345,
    ...) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(sce)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'sce'")
        }
    } else {
        sample <- rep(1, ncol(sce))
    }

    message(paste0(date(), " ... Running 'cxds_bcds_hybrid'"))

    ## Define result matrix for all samples
    if ("estNdbl" %in% names(list(...))) {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            hybrid_score = numeric(ncol(sce)),
            hybrid_call = logical(ncol(sce)))
    } else {
        output <- S4Vectors::DataFrame(row.names = colnames(sce),
            hybrid_score = numeric(ncol(sce)))
    }

    ## Loop through each sample and run cxds_bcds_hybrid
    samples <- unique(sample)
    for (i in seq_len(length(samples))) {
        sceSampleInd <- sample == samples[i]
        sceSample <- sce[, sceSampleInd]

        counts(sceSample) <- as(counts(sceSample), "dgCMatrix")
        result <- withr::with_seed(seed, scds::cxds_bcds_hybrid(sce = sceSample,
            ...))

        if ("hybrid_call" %in% colnames(SummarizedExperiment::colData(result))) {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("hybrid_score", "hybrid_call")]
        } else {
            output[sceSampleInd, ] <- SummarizedExperiment::colData(result)[,
                c("hybrid_score")]
        }
    }

    colnames(output) <- paste0("scds_", colnames(output))
    colData(sce) = cbind(colData(sce), output)

    return(sce)
}
