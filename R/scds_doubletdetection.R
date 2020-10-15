
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
#' @param ntop See \link[scds]{cxds} for more information. Default \code{500}.
#' @param binThresh See \link[scds]{cxds} for more information. Default \code{0}.
#' @param verb See \link[scds]{cxds} for more information. Default \code{FALSE}.
#' @param retRes See \link[scds]{cxds} for more information. Default \code{FALSE}.
#' @param estNdbl See \link[scds]{cxds} for more information. Default \code{FALSE}.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{cxds_score} and optionally \emph{cxds_call}.
#'  Please refer to the documentation of \link[scds]{cxds} for details.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runCxds(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment counts counts<-
runCxds <- function(inSCE,
    sample = NULL,
    seed = 12345,
    ntop = 500,
    binThresh = 0,
    verb = FALSE,
    retRes = FALSE,
    estNdbl = FALSE,
    useAssay = "counts") {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'cxds'"))

    ## Getting current arguments
    #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
    argsList <- mget(names(formals()),sys.frame(sys.nframe()))

    ## Define result matrix for all samples
    if (estNdbl) {
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

        mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
        counts(sceSample) <- .convertToMatrix(mat)

        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed, scds::cxds(sce = sceSample,
                                                           ntop = nGene,
                                                           binThresh = binThresh,
                                                           verb = verb,
                                                           retRes = retRes,
                                                           estNdbl = estNdbl))}, silent = TRUE)
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

    inSCE@metadata$runCxds <- argsList[-1]
    inSCE@metadata$runCxds$packageVersion <- utils::packageDescription("scds")$Version

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
#' @param ntop See \link[scds]{bcds} for more information. Default \code{500}.
#' @param srat See \link[scds]{bcds} for more information. Default \code{1}.
#' @param verb See \link[scds]{bcds} for more information. Default \code{FALSE}.
#' @param retRes See \link[scds]{bcds} for more information. Default \code{FALSE}.
#' @param nmax See \link[scds]{bcds} for more information. Default \code{"tune"}.
#' @param varImp See \link[scds]{bcds} for more information. Default \code{FALSE}.
#' @param estNdbl See \link[scds]{bcds} for more information. Default \code{FALSE}.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{bcds} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{bcds_score} and optionally \emph{bcds_call}.
#'  Please refer to the documentation of \link[scds]{bcds} for details.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runBcds(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment counts counts<-
runBcds <- function(inSCE,
    sample = NULL,
    seed = 12345,
    ntop = 500,
    srat = 1,
    verb = FALSE,
    retRes = FALSE,
    nmax = "tune",
    varImp = FALSE,
    estNdbl = FALSE,
    useAssay = "counts"
    ) {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'bcds'"))

    ## Getting current arguments
    #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
    argsList <- mget(names(formals()),sys.frame(sys.nframe()))

    ## Define result matrix for all samples
    if (estNdbl) {
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

        mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
        counts(sceSample) <- .convertToMatrix(mat)

        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed,
            scds::bcds(sce = sceSample,
                       ntop = nGene,
                       srat = srat,
                       verb = verb,
                       retRes = retRes,
                       nmax = nmax,
                       varImp = varImp,
                       estNdbl = estNdbl
            ))}, silent = TRUE)
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

    inSCE@metadata$runBcds <- argsList[-1]
    inSCE@metadata$runBcds$packageVersion <- utils::packageDescription("scds")$Version

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
#' @param nTop The number of top varialbe genes to consider. Used in both \code{csds}
#' and \code{bcds}. Default \code{500}.
#' @param cxdsArgs See \link[scds]{cxds_bcds_hybrid} for more information. Default \code{NULL}.
#' @param bcdsArgs See \link[scds]{cxds_bcds_hybrid} for more information. Default \code{NULL}.
#' @param verb See \link[scds]{cxds_bcds_hybrid} for more information. Default \code{FALSE}.
#' @param estNdbl See \link[scds]{cxds_bcds_hybrid} for more information. Default \code{FALSE}.
#' @param force See \link[scds]{cxds_bcds_hybrid} for more information. Default \code{FALSE}.
#' @param useAssay  A string specifying which assay in the SCE to use.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  \link[scds]{cxds_bcds_hybrid} output appended to the
#'  \link[SummarizedExperiment]{colData} slot. The columns include
#'  \emph{hybrid_score} and optionally \emph{hybrid_call}.
#'  Please refer to the documentation of \link[scds]{cxds_bcds_hybrid} for
#'  details.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runCxdsBcdsHybrid(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment counts counts<-
runCxdsBcdsHybrid <- function(inSCE,
    sample = NULL,
    seed = 12345,
    nTop = 500,
    cxdsArgs = list(),
    bcdsArgs = list(),
    verb = FALSE,
    estNdbl = FALSE,
    force = FALSE,
    useAssay = "counts") {

    if (!is.null(sample)) {
        if (length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number",
                " of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    message(paste0(date(), " ... Running 'cxds_bcds_hybrid'"))

    ## Getting current arguments
    #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
    argsList <- mget(names(formals()),sys.frame(sys.nframe()))

    ## Define result matrix for all samples
    if (estNdbl) {
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

        mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
        counts(sceSample) <- .convertToMatrix(mat)

        result <- NULL
        nGene <- 500
        while(!inherits(result, "SingleCellExperiment") & nGene > 0) {
          try({result <- withr::with_seed(seed, scds::cxds_bcds_hybrid(sce = sceSample,
                                                                       cxdsArgs=c(list(ntop = nGene), cxdsArgs),
                                                                       bcdsArgs=c(list(ntop = nGene), bcdsArgs),
                                                                       verb = verb,
                                                                       estNdbl = estNdbl,
                                                                       force = force))}, silent = TRUE)
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

    inSCE@metadata$runCxdsBcdsHybrid <- argsList[-1]
    inSCE@metadata$runCxdsBcdsHybrid$packageVersion <- utils::packageDescription("scds")$Version

    return(inSCE)
}
