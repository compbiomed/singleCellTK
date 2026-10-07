#' @title Calculate per-cell QC metrics
#' @description Calculates general quality control metrics for each cell in
#' the count matrix: total counts, detected features, percent of counts in
#' the top features, and the same metrics for gene sets (such as mitochondrial
#' genes) and alternative experiments. The metrics and column names are the
#' same as those of \code{scater::addPerCellQC}, which is deprecated.
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param useAssay A string specifying which assay in the SCE to use. Default
#' \code{"counts"}.
#' @param mitoRef Character. The species used to extract mitochondrial genes ID 
#' from build-in mitochondrial geneset in SCTK. Available species options are 
#' \code{"human"} and \code{"mouse"}. Default is \code{"human"}. 
#' @param mitoIDType Character. Types of mitochondrial gene id. SCTK supports 
#' \code{"symbol"}, \code{"entrez"}, \code{"ensembl"} and 
#' \code{"ensemblTranscriptID"}. It is used with \code{mitoRef} to extract 
#' mitochondrial genes from build-in mitochondrial geneset in SCTK. Default 
#' \code{NULL}. 
#' @param mitoGeneLocation Character. Describes the location within \code{inSCE}
#' where the gene identifiers in the mitochondrial gene sets should be located.
#' If set to \code{"rownames"} then the features will be searched for among 
#' \code{rownames(inSCE)}. This can also be set to one of the column names of 
#' \code{rowData(inSCE)} in which case the gene identifies will be mapped to 
#' that column in the \code{rowData} of \code{inSCE}. See 
#' \code{\link{featureIndex}} for more information. If this parameter is set to
#' \code{NULL}, then no mitochondrial metrics will be calculated.
#' Default \code{"rownames"}.
#' @param mitoPrefix Character. The prefix used to get mitochondrial gene from 
#' either \code{rownames(inSCE)} or columns of \code{rowData(inSCE)} specified 
#' by \code{mitoGeneLocation}. This parameter is usually used to extract mitochondrial 
#' genes from the gene symbol. For example, \code{mitoPrefix = "^MT-"} can be used 
#' to detect mito gene symbols like "MT-ND4". Note that case is ignored so "mt-"
#' will still match "MT-ND4". Default \code{"^MT-"}.
#' @param mitoID Character. A vector of mitochondrial genes to be quantified.  
#' @param collectionName Character. Name of a \code{GeneSetCollection} obtained
#' by using one of the \code{importGeneSet*} functions. Default \code{NULL}.
#' @param geneSetList List of gene sets to be quantified. The genes in the
#' assays will be matched to the genes in the list based on
#' \code{geneSetListLocation}. Default \code{NULL}.
#' @param geneSetListLocation Character or numeric vector. If set to 
#' \code{'rownames'}, then the genes in \code{geneSetList} will be looked up in 
#' \code{rownames(inSCE)}. If another character is supplied, then genes will be 
#' looked up in the column names of \code{rowData(inSCE)}. A character vector 
#' with the same length as \code{geneSetList} can be supplied if the IDs for 
#' different gene sets are found in different places, including a mixture of 
#' \code{'rownames'} and \code{rowData(inSCE)}. An integer or integer vector can
#' be supplied to denote the column index in \code{rowData(inSCE)}. Default 
#' \code{'rownames'}.
#' @param geneSetCollection Class of \code{GeneSetCollection} from package
#' GSEABase. The location of the gene IDs in \code{inSCE} should be in the
#' \code{description} slot of each gene set and should follow the
#' same notation as \code{geneSetListLocation}. The function 
#' \code{\link[GSEABase]{getGmt}} can be used to read in gene sets from a GMT 
#' file. If reading a GMT file, the second column for each gene set should be 
#' the description denoting the location of the gene IDs in \code{inSCE}. These 
#' gene sets will be included with those from \code{geneSetList} if both 
#' parameters are provided.
#' @param percent_top An integer vector. Each element is treated as a number of
#' top genes to compute the percentage of library size occupied by the most 
#' highly expressed genes in each cell. Default \code{c(50, 100, 200, 500)}.
#' @param use_altexps Logical scalar indicating whether QC statistics should
#' be computed for alternative Experiments in \code{inSCE} 
#' (\code{altExps(inSCE)}). If \code{TRUE}, statistics are computed for all 
#' alternative experiments. Alternatively, an integer or character vector 
#' specifying the alternative Experiments to use to compute QC statistics. 
#' Alternatively \code{NULL}, in which case alternative experiments are not 
#' used. Default \code{FALSE}.
#' @param flatten Logical scalar indicating whether the nested
#' \link[S4Vectors]{DataFrame-class} in the output should be flattened. Default 
#' \code{TRUE}.
#' @param detectionLimit A numeric scalar specifying the lower detection limit
#' for expression. Default \code{0}
#' @param BPPARAM A \link{BiocParallelParam} object specifying whether the QC
#' calculations should be parallelized. Default 
#' \code{BiocParallel::SerialParam()}.
#' @details 
#' This function allows multiple ways to import mitochondrial genes and quantify 
#' their expression in cells. \code{mitoGeneLocation} is required for all
#' methods to point to the location within inSCE object that stores the
#' mitochondrial gene IDs or Symbols. The various ways mito genes can be 
#' specified are:

#' \itemize{
#'   \item A combination of \code{mitoRef} and \code{mitoIDType} 
#'   parameters can be used to load pre-built mitochondrial gene sets stored
#'   in the SCTK package. These parameters are used in the
#'   \link{importMitoGeneSet} function.
#'   \item The \code{mitoPrefix} parameter can be used to search for features
#'   matching a particular pattern. The default pattern is an "MT-" 
#'   at the beginning of the ID.
#'   \item The \code{mitoID} parameter can be used to directy supply a vector of 
#'   mitochondrial gene IDs or names. Only features that exactly match items
#'   in this vector will be included in the mitochondrial gene set. 
#' }
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#' cell QC metrics added to the \link{colData} slot. 
#' @seealso
#' \code{link{plotRunPerCellQCResults}}, \code{\link{runCellQC}}
#' @examples
#' data(scExample, package = "singleCellTK")
#' mito.ix = grep("^MT-", rowData(sce)$feature_name)
#' geneSet <- list("Mito"=rownames(sce)[mito.ix])
#' sce <- runPerCellQC(sce, geneSetList = geneSet)
#' @export
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom S4Vectors metadata metadata<-
runPerCellQC <- function(inSCE,
                         useAssay = "counts",
                         mitoGeneLocation = "rownames",                         
                         mitoRef = c(NULL, "human", "mouse"),
                         mitoIDType = c("ensembl", "symbol", "entrez", "ensemblTranscriptID"),
                         mitoPrefix = "MT-",
                         mitoID = NULL,
                         collectionName = NULL,
                         geneSetList = NULL,
                         geneSetListLocation = "rownames",
                         geneSetCollection = NULL,
                         percent_top = c(50, 100, 200, 500),
                         use_altexps = FALSE,
                         flatten = TRUE,
                         detectionLimit = 0,
                         BPPARAM = BiocParallel::SerialParam()
) {

  p <- paste0(date(), " ... Running 'perCellQCMetrics'")
  message(p)
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))
  mitoRef <- match.arg(mitoRef)
  mitoIDType <- match.arg(mitoIDType)
  
  ## Add mito gene collection from built-in mito gene sets in SCTK package
  if (!is.null(mitoGeneLocation) && mitoGeneLocation == "rownames") {
    features <- rownames(inSCE)
  } else if (mitoGeneLocation %in% names(rowData(inSCE))) {
    features <- rowData(inSCE)[[mitoGeneLocation]]
  } else {
    stop("'mitoGeneLocation' needs to be set to 'rownames' or a column name of 'rowData(inSCE)'")
  }
  
  if(!is.null(mitoGeneLocation)) {
    mitoGS <- NULL
    
    # Getting rid of any previous gene sets with the same name
    metadata(inSCE)$sctk$genesets$mito <- NULL

    ## Add mito genes specified in mitoID
    if (!is.null(mitoID)) {
      message(date(), " ...... Attempting to find mitochondiral by identifing features in '", mitoGeneLocation, "' that match those given in `mitoID`.")
      inSCE <- importGeneSetsFromList(inSCE, geneSetList = list(mito=mitoID), 
                                      collectionName = "mito", 
                                      by = mitoGeneLocation,
                                      noMatchError = FALSE)
      mitoGS <- metadata(inSCE)$sctk$genesets$mito[[1]]
    }
    
    if (!is.null(mitoRef) & !is.null(mitoIDType) & is.null(mitoGS)) {
      message(date(), " ...... Attempting to find mitochondrial genes by identifying features in '", mitoGeneLocation, 
      "' that match mitochondrial genes from reference '", mitoRef, "' and ID type '", mitoIDType, "'.")
      
      inSCE <- importMitoGeneSet(inSCE, reference = mitoRef, id = mitoIDType, 
                                 by = mitoGeneLocation, collectionName = "mito",
                                 noMatchError = FALSE)
      mitoGS <- metadata(inSCE)$sctk$genesets$mito[[1]]
    }

    ## Add mito genes by greping mitoPrefix from the mitoGeneLocation
    if (!is.null(mitoPrefix) & is.null(mitoGS)) {
      message(date(), " ...... Attempting to find mitochondrial genes by identifying features in '", mitoGeneLocation, "' starting with '", mitoPrefix, "'")
      mitoG <- grep(paste0("^", mitoPrefix), features, value = TRUE, ignore.case = TRUE)
      
      if(length(mitoG) > 0) {
        inSCE <- importGeneSetsFromList(inSCE, list(mito=mitoG), 
                                        by = mitoGeneLocation, collectionName = "mito",
                                        noMatchError = FALSE)
        mitoGS <- metadata(inSCE)$sctk$genesets$mito[[1]]
      }
    }

    # Creating a new gene set collection from mitoGS and any other user supplied geneSetCollections
    if (is.null(mitoGS)) {
      message(date(), " ...... No mitochondrial genes found. Skipping quantification of mitochondrial metrics.")      
    } else {
      if (!is.null(geneSetCollection)) {
        geneSetCollection <- GSEABase::GeneSetCollection(
          c(mitoGS, unlist(geneSetCollection))
        )
      } else {
        geneSetCollection <- GSEABase::GeneSetCollection(mitoGS)
      }    
    }
  }

  ## Add GeneSetColletion that has been previously imported
  if(!is.null(collectionName)) {
    gsc <- .retrieveGeneSetCollection(inSCE = inSCE,
                                      collectionName = collectionName)
    geneSetList <- c(geneSetList, GSEABase::geneIds(gsc))
  }

  ## Add gene sets in 'geneSetCollection' to 'geneSetList', if available
  original.length <- length(geneSetList)
  geneSetCollectionLocation <- c()
  if(!is.null(geneSetCollection)) {

    ## Get the location where the gene set Ids are stored in SCE object
    geneSetCollectionLocation <- vapply(geneSetCollection, 
                                        GSEABase::description,
                                        FUN.VALUE = character(1))
    #character(length(names(geneSetCollection)))

    ## If blank/null/NA, then set to rownames by default
    ix <- geneSetCollectionLocation == "" ||
      is.na(geneSetCollectionLocation) ||
      is.null(geneSetCollectionLocation)
    geneSetCollectionLocation[ix] <- "rownames"

    ## Add gene sets to geneSetList
    geneSetList <- c(geneSetList, GSEABase::geneIds(geneSetCollection))
  }

  ## Set up locations for all lists
  if(length(geneSetListLocation) == 1) {
    locations <- c(rep(geneSetListLocation, original.length),
                   geneSetCollectionLocation)
  } else if(length(geneSetListLocation) == original.length) {
    locations <- c(geneSetListLocation, geneSetCollectionLocation)
  } else {
    stop("'geneSetListLocation' needs to be of length 1 and will be used for ",
         "all gene sets in 'geneSetList' or needs to be the length of ",
         "'geneSetList'.")
  }

  ## Process gene sets in the list
  geneSets <- list()
  if(!is.null(geneSetList)) {

    if(is.null(names(geneSetList))) {
      names(geneSetList) <- paste0("GeneSet_", seq_along(geneSetList))
    }

    for(i in seq_along(geneSetList)) {
      if(locations[i] == "rownames") {
        ## IDs located in rownames
        gs.i <- intersect(geneSetList[[i]], rownames(inSCE))
        gs.diff <- setdiff(geneSetList[[i]], rownames(inSCE))

        if(length(gs.i) > 0) {
          gs.index <- which(rownames(inSCE) %in% gs.i)
          geneSets[[names(geneSetList)[i]]] <- rownames(inSCE)[gs.index]
        } else {
          warning("No features for '", names(geneSetList)[i], 
                  "' were found in 'rownames(inSCE)'. Excluding this gene set.")
        }

        if(length(gs.diff) > 0 & length(gs.i) > 0) {
          warning("Some features were not found in 'rownames(inSCE)' for ", 
                  "gene set '", names(geneSetList)[i], "': ",
                  paste(gs.diff, collapse=","))
        }

      } else {
        ## If IDs located in rowData
        ## Check for existence of location in rowData(inSCE)
        if(any(locations[i] == seq(ncol(rowData(inSCE))))) {
          temp.location <- as.numeric(locations[i])
        } else if (any(locations[i] == colnames(rowData(inSCE)))) {
          temp.location <- as.character(locations[i])
        } else {
          stop("'", locations[i], "' was not found in the column names of ", 
               "'rowData(inSCE).'")
        }

        gs.i <- intersect(geneSetList[[i]], rowData(inSCE)[,temp.location])
        gs.diff <- setdiff(geneSetList[[i]], rowData(inSCE)[,temp.location])

        if(length(gs.i) > 0) {
          gs.index <- which(rowData(inSCE)[,temp.location] %in% gs.i)
          geneSets[[names(geneSetList)[i]]] <- rownames(inSCE)[gs.index]
        } else {
          p <- paste0("No features for gene set '", names(geneSetList)[i],
                      "' were found in 'rowData(inSCE)' under column '",
                      temp.location, "'. Excluding this gene set.", sep="")
          warning(p)
        }

        if(length(gs.diff) > 0 & length(gs.i) > 0) {
          warning(paste0("Some features were not found in 'rowData(inSCE)' ", 
                         "under column '", temp.location, " 'for gene set '",
                         names(geneSetList)[i], "': ",
                         paste(gs.diff, collapse=","), sep=""))
        }
      }
    }
  }

  if(length(geneSets) == 0) {
    geneSets <- NULL
  }
  
  colData(inSCE)$sum <- NULL
  colData(inSCE)$detected <- NULL
  colData(inSCE)$percent.top_50 <- NULL
  colData(inSCE)$percent.top_100 <- NULL
  colData(inSCE)$percent.top_200 <- NULL
  colData(inSCE)$percent.top_500 <- NULL
  colData(inSCE)$total <- NULL
    
  qcMetrics <- .perCellQCMetrics(inSCE, useAssay = useAssay,
                                 subsets = geneSets,
                                 percentTop = percent_top,
                                 useAltExps = use_altexps,
                                 flatten = flatten,
                                 detectionLimit = detectionLimit)
  colData(inSCE) <- cbind(colData(inSCE), qcMetrics)

  ## rename mito gene columns in colData(inSCE)
  names(colData(inSCE)) <- gsub('subsets_mito', 'mito', names(colData(inSCE)))

  argsList <- argsList[!names(argsList) %in% ("BPPARAM")]
  metadata(inSCE)$sctk$runPerCellQC$all_cells <- argsList[-1]
  metadata(inSCE)$sctk$runPerCellQC$all_cells$packageVersion <- 
    utils::packageDescription("singleCellTK")$Version

  if(is.null(geneSets)){
    geneSets <- as.character(geneSets)
  }
  metadata(inSCE)$sctk$runPerCellQC$all_cells$geneSets <- geneSets

  return(inSCE)
}

# Per-cell QC metrics with the same columns and values as
# scuttle::perCellQCMetrics(), which is deprecated: total counts, detected
# features, percent of counts in the top N features, the same three metrics
# for each subset of features and each alternative experiment, and the total
# over the main and alternative experiments. Returns a DataFrame with one row
# per cell, nested like scuttle's when flatten = FALSE.
.perCellQCMetrics <- function(inSCE, useAssay = "counts", subsets = NULL,
                              percentTop = integer(0), useAltExps = FALSE,
                              flatten = TRUE, detectionLimit = 0,
                              chunkSize = 5000) {
  percentTop <- sort(as.integer(percentTop))
  altNames <- .selectAltExps(inSCE, useAltExps, useAssay)
  if (is.null(subsets)) subsets <- list()

  main <- .countMetrics(assay(inSCE, useAssay), detectionLimit, percentTop,
                        subsets, chunkSize)
  libSize <- main$sum
  altMetrics <- lapply(altNames, function(name) {
    alt <- SingleCellExperiment::altExp(inSCE, name)
    .countMetrics(assay(alt, useAssay), detectionLimit, integer(0), list(),
                  chunkSize)
  })
  total <- libSize + Reduce(`+`, lapply(altMetrics, `[[`, "sum"),
                            numeric(length(libSize)))

  cells <- colnames(inSCE)
  topPercent <- main$top / libSize * 100
  colnames(topPercent) <- as.character(percentTop)
  subsetDF <- S4Vectors::make_zero_col_DFrame(ncol(inSCE))
  for (name in names(subsets)) {
    subsetDF[[name]] <- S4Vectors::DataFrame(
      sum = main$subsets[[name]]$sum,
      detected = main$subsets[[name]]$detected,
      percent = main$subsets[[name]]$sum / libSize * 100
    )
  }
  altDF <- S4Vectors::make_zero_col_DFrame(ncol(inSCE))
  for (i in seq_along(altNames)) {
    altDF[[altNames[i]]] <- S4Vectors::DataFrame(
      sum = altMetrics[[i]]$sum,
      detected = altMetrics[[i]]$detected,
      percent = altMetrics[[i]]$sum / total * 100
    )
  }
  out <- S4Vectors::DataFrame(sum = libSize, detected = main$detected,
                              row.names = cells)
  out$percent.top <- topPercent
  out$subsets <- subsetDF
  out$altexps <- altDF
  out$total <- total
  if (isTRUE(flatten)) out <- .flattenQCMetrics(out)
  out
}

# Names of the alternative experiments to include in QC metrics. useAltExps
# is TRUE (all), FALSE (none), NULL (those with an assay named useAssay,
# scuttle's default), or names or indices.
.selectAltExps <- function(inSCE, useAltExps, useAssay) {
  allAlt <- SingleCellExperiment::altExpNames(inSCE)
  if (isTRUE(useAltExps)) return(allAlt)
  if (isFALSE(useAltExps)) return(character(0))
  if (!is.null(useAltExps)) return(allAlt[useAltExps])
  hasAssay <- vapply(allAlt, function(n) {
    alt <- SingleCellExperiment::altExp(inSCE, n)
    useAssay %in% SummarizedExperiment::assayNames(alt)
  }, logical(1))
  allAlt[hasAssay]
}

# Column sums, detected counts, top-N sums, and per-subset sums and detected
# counts of a features x cells matrix, processed in column chunks so that
# sparse and DelayedArray inputs are never fully densified.
.countMetrics <- function(mat, detectionLimit, percentTop, subsets,
                          chunkSize) {
  nCells <- ncol(mat)
  sums <- numeric(nCells)
  detected <- numeric(nCells)
  top <- matrix(0, nrow = nCells, ncol = length(percentTop))
  subsetIndex <- lapply(subsets, function(s) {
    idx <- seq_len(nrow(mat))
    names(idx) <- rownames(mat)
    unname(idx[s])
  })
  subsetSums <- lapply(subsets, function(s) numeric(nCells))
  subsetDetected <- lapply(subsets, function(s) numeric(nCells))
  for (start in seq(1, nCells, by = chunkSize)) {
    cols <- seq(start, min(start + chunkSize - 1, nCells))
    chunk <- methods::as(mat[, cols, drop = FALSE], "CsparseMatrix")
    chunk <- methods::as(chunk, "dMatrix")
    sums[cols] <- Matrix::colSums(chunk)
    detected[cols] <- Matrix::colSums(chunk > detectionLimit)
    if (length(percentTop) > 0) {
      top[cols, ] <- .topSums(chunk, percentTop)
    }
    for (name in names(subsets)) {
      part <- chunk[subsetIndex[[name]], , drop = FALSE]
      subsetSums[[name]][cols] <- Matrix::colSums(part)
      subsetDetected[[name]][cols] <- Matrix::colSums(part > detectionLimit)
    }
  }
  subsetMetrics <- lapply(stats::setNames(nm = names(subsets)), function(n) {
    list(sum = subsetSums[[n]], detected = subsetDetected[[n]])
  })
  list(sum = sums, detected = detected, top = top, subsets = subsetMetrics)
}

# Sum of the largest N values in each column of a sparse matrix, for each N
# in percentTop. Values that are not stored are zero, so only stored values
# need sorting unless the matrix has negative values.
.topSums <- function(chunk, percentTop) {
  if (any(chunk@x < 0)) {
    dense <- as.matrix(chunk)
    return(t(apply(dense, 2, function(v) {
      v <- sort(v, decreasing = TRUE)
      vapply(percentTop, function(n) sum(v[seq_len(min(n, length(v)))]),
             numeric(1))
    })))
  }
  stored <- diff(chunk@p)
  column <- rep.int(seq_len(ncol(chunk)), stored)
  ord <- order(column, -chunk@x)
  values <- chunk@x[ord]
  rank <- sequence(stored)
  vapply(percentTop, function(n) {
    keep <- rank <= n
    out <- numeric(ncol(chunk))
    s <- rowsum(values[keep], column[ord][keep])
    out[as.integer(rownames(s))] <- s[, 1]
    out
  }, numeric(ncol(chunk)))
}

# Flatten nested QC metrics into scuttle's column names, for example
# percent.top_50, subsets_mito_sum, and altexps_ERCC_percent.
.flattenQCMetrics <- function(qc) {
  out <- qc[, c("sum", "detected"), drop = FALSE]
  top <- qc$percent.top
  for (n in colnames(top)) out[[paste0("percent.top_", n)]] <- top[, n]
  for (group in c("subsets", "altexps")) {
    nested <- qc[[group]]
    for (name in colnames(nested)) {
      for (metric in colnames(nested[[name]])) {
        out[[paste(group, name, metric, sep = "_")]] <- nested[[name]][[metric]]
      }
    }
  }
  out$total <- qc$total
  out
}
