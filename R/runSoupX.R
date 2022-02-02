#' @title Detecting contamination with SoupX
#' @description A wrapper function for \link[SoupX]{adjustCounts}. Identify
#' potential contamination from experimental factors such as ambient RNA. Visit
#' \href{https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html}{their vignette}
#' for better understanding. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object.
#' @param sample A single character specifying a name that can be found in
#' \code{colData(inSCE)} to directly use the cell annotation; or a character
#' vector with as many elements as cells to indicates which sample each cell
#' belongs to. SoupX will be run on cells from each sample separately. Default 
#' \code{NULL}. 
#' @param useAssay A single character string specifying which assay in 
#' \code{inSCE} to use. Default \code{'counts'}.
#' @param background A numeric matrix of counts or a 
#' \linkS4class{SingleCellExperiment} object with the matrix in \code{assay} 
#' slot. It should have the same structure as \code{inSCE} except it contains 
#' the matrix including empty droplets. Default \code{NULL}.
#' @param bgAssayName A single character string specifying which assay in 
#' \code{background} to use when \code{background} is a 
#' \linkS4class{SingleCellExperiment} object. If \code{NULL}, the function
#' will use the same value as \code{useAssay}. Default \code{NULL}.
#' @param bgBatch The same thing as \code{sample} but for \code{background}. Can
#' be a single character only when \code{background} is a 
#' \linkS4class{SingleCellExperiment} object. Default \code{NULL}.
#' @param assayName A single character string of the output corrected matrix. 
#' Default \code{"SoupX"} when not using a background, otherwise, 
#' \code{"SoupX_bg"}.
#' @param cluster Prior knowledge of clustering labels on cells. A single 
#' character string for specifying clustering label stored in 
#' \code{colData(inSCE)}, or a character vector with as many elements as cells.
#' When not supplied, SCTK default clustering method will be applied (see 
#' \href{https://camplab.net/sctk/v2.4.1/articles/articles/console_analysis_tutorial.html#normalization}{detail}). 
#' @param tfidfMin SoupX parameter for \link[SoupX]{autoEstCont}. Minimum value 
#' of tfidf to accept for a marker gene. Default \code{1}. See 
#' \code{?SoupX::autoEstCont}. 
#' @param soupQuantile SoupX parameter for \link[SoupX]{autoEstCont}. Only use 
#' genes that are at or above this expression quantile in the soup. This 
#' prevents inaccurate estimates due to using genes with poorly constrained 
#' contribution to the background. Default \code{0.9}. See 
#' \code{?SoupX::autoEstCont}. 
#' @param adjustMethod SoupX parameter for \link[SoupX]{adjustCounts}. Method to
#' use for correction. One of \code{'multinomial'}, \code{'soupOnly'}, or 
#' \code{'subtraction'}. See \code{?SoupX::adjustCounts}. 
#' @param roundToInt SoupX parameter for \link[SoupX]{adjustCounts}. Should the 
#' resulting matrix be rounded to integers? Default \code{FALSE}. See 
#' \code{?SoupX::adjustCounts}. 
#' @param pCut SoupX parameter for \link[SoupX]{adjustCounts}. The p-value 
#' cut-off used when \code{method = 'soupOnly'}. Default \code{0.01}. See 
#' \code{?SoupX::adjustCounts}. 
#' @return The input \code{inSCE} object with \code{soupX_nUMIs}, 
#' \code{soupX_clustrers}, \code{soupX_contamination} added to \code{colData} 
#' slot; \code{soupX_{sample}_est} and \code{soupX_{sample}_counts} for each 
#' sample added in \code{rowData} slot; and other computation metrics at 
#' \code{getSoupX(inSCE)}. Replace "soupX" to "soupX_bg" when \code{background}
#' is used. 
#' @export
#' @author Yichen Wang
runSoupX <- function(inSCE, 
                     sample = NULL,
                     useAssay = "counts", 
                     background = NULL, 
                     bgAssayName = NULL,
                     bgBatch = NULL,
                     assayName = ifelse(is.null(background), 
                                        "SoupX", "SoupX_bg"), 
                     cluster = NULL, 
                     tfidfMin = 1, 
                     soupQuantile = 0.9,
                     adjustMethod = c("subtraction", "soupOnly", "multinomial"),
                     roundToInt = FALSE,
                     pCut = 0.01) {
    SingleBGBatchForAllBatch <- FALSE
    adjustMethod <- match.arg(adjustMethod)
    if(!is.null(sample)) {
        if (length(sample) == 1) {
            if (!sample %in% names(SummarizedExperiment::colData(inSCE))) {
                stop("Specified Sample variable not found in colData.")
            }
            sample <- SummarizedExperiment::colData(inSCE)[[sample]]
        } else if(length(sample) != ncol(inSCE)) {
            stop("'sample' must be the same length as the number ", 
                 "of columns in 'inSCE'")
        }
        if (is.factor(sample)) {
            sample <- as.character(sample)
        }
        uniqSample <- unique(sample)
        if (!is.null(background)) {
            if (is.null(bgBatch)) {
                # 'sample' specified, 'bgBatch' not
                bgBatch <- rep("all_cells", ncol(background))
                SingleBGBatchForAllBatch <- TRUE
            } else {
                # 'sample' & 'bgBatch' both specified
                if (length(bgBatch) == 1) {
                    if (!bgBatch %in% names(SummarizedExperiment::colData(background))) {
                        stop("Specified bgBatch variable not found ", 
                             "in background colData")
                    }
                    bgBatch <- SummarizedExperiment::colData(background)[[bgBatch]]
                } else if(length(bgBatch) != ncol(background)) {
                    stop("'bgBatch' must be the same length as ",
                         "the number of columns in 'background'")
                }
                if (is.factor(bgBatch)) {
                    bgBatch <- as.character(bgBatch)
                }
                # Check 'sample's in cell matrix can all be found in bgBatch
                if (!all(sample) %in% bgBatch) {
                    stop("Not all samples can be found in 'bgBatch'.")
                }
            }
        }
    } else {
        # 'sample' not specified. 
        sample <- rep("all_cells", ncol(inSCE))
        if (!is.null(background)) {
            if (!is.null(bgBatch)) {
                warning("Using all background because 'sample' not specified.")
            }
            bgBatch <- rep("all_cells", ncol(background))
        }
        uniqSample <- "all_cells"
    }
    
    message(paste0(date(), " ... Running 'SoupX'"))
    
    results <- list()
    sampleIdx <- list()
    for (s in uniqSample) {
        message(paste0(date(), " ... Running 'SoupX' on sample: ", s))
        cellIdx <- sample == s
        sampleIdx[[s]] <- cellIdx
        tempSCE <- inSCE[,cellIdx]
        if (isTRUE(SingleBGBatchForAllBatch)) {
            res <- .SoupXOneBatch(inSCE = tempSCE, 
                                  useAssay = useAssay,
                                  background = background,
                                  bgAssayName = bgAssayName,
                                  cluster = cluster,
                                  tfidfMin = tfidfMin,
                                  soupQuantile = soupQuantile,
                                  adjustMethod = adjustMethod,
                                  roundToInt = roundToInt,
                                  pCut = pCut)
        } else {
            if (!is.null(background)) {
                bgIdx <- bgBatch == s
                tempBG <- background[,bgIdx]
            } else {
                tempBG <- NULL
            }
            res <- .SoupXOneBatch(inSCE = tempSCE, 
                                  useAssay = useAssay,
                                  background = tempBG,
                                  bgAssayName = bgAssayName,
                                  cluster = cluster,
                                  tfidfMin = tfidfMin,
                                  soupQuantile = soupQuantile,
                                  adjustMethod = adjustMethod,
                                  roundToInt = roundToInt,
                                  pCut = pCut)
        }
        results[[s]] <- res
    }
    # Initiate new assay by copying the input selection
    # And then replace with new value at sample indices
    # Similarly for colData and rowData
    corrAssay <- SummarizedExperiment::assay(inSCE, useAssay)
    newColData <- SummarizedExperiment::colData(inSCE)
    if (!is.null(background)) {
        newColData$soupX_bg_nUMIs <- NA
        newColData$soupX_bg_clusters <- NA
        newColData$soupX_bg_contamination <- NA
    } else {
        newColData$soupX_nUMIs <- NA
        newColData$soupX_clusters <- NA
        newColData$soupX_contamination <- NA
    }
    newRowData <- SummarizedExperiment::rowData(inSCE)
    for (s in names(results)) {
        corrAssay[,sampleIdx[[s]]] <- results[[s]]$out
        if (!is.null(background)) {
            newColData$soupX_bg_nUMIs[sampleIdx[[s]]] <- results[[s]]$sc$metaData$nUMIs
            newColData$soupX_bg_clusters[sampleIdx[[s]]] <- results[[s]]$sc$metaData$clusters
            newColData$soupX_bg_contamination[sampleIdx[[s]]] <- results[[s]]$sc$metaData$rho
            newRowData[[paste0("soupX_bg_",s,"_est")]] <- results[[s]]$sc$soupProfile$est
            newRowData[[paste0("soupX_bg_",s,"_counts")]] <- results[[s]]$sc$soupProfile$counts
            S4Vectors::metadata(inSCE)$sctk$SoupX_bg[[s]] <- results[[s]]$sc$fit
            S4Vectors::metadata(inSCE)$sctk$SoupX_bg[[s]]$nDropUMIs <- results[[s]]$sc$nDropUMIs
        } else {
            newColData$soupX_nUMIs[sampleIdx[[s]]] <- results[[s]]$sc$metaData$nUMIs
            newColData$soupX_clusters[sampleIdx[[s]]] <- results[[s]]$sc$metaData$clusters
            newColData$soupX_contamination[sampleIdx[[s]]] <- results[[s]]$sc$metaData$rho
            newRowData[[paste0("soupX_",s,"_est")]] <- results[[s]]$sc$soupProfile$est
            newRowData[[paste0("soupX_",s,"_counts")]] <- results[[s]]$sc$soupProfile$counts
            S4Vectors::metadata(inSCE)$sctk$SoupX[[s]] <- results[[s]]$sc$fit
            S4Vectors::metadata(inSCE)$sctk$SoupX[[s]]$nDropUMIs <- results[[s]]$sc$nDropUMIs
        }
    }
    expData(inSCE, assayName, tag = "raw") <- corrAssay
    SummarizedExperiment::colData(inSCE) <- newColData
    SummarizedExperiment::rowData(inSCE) <- newRowData
    
    return(inSCE)
}

.SoupXOneBatch <- function(inSCE, 
                           useAssay, 
                           background, 
                           bgAssayName,
                           cluster, 
                           tfidfMin, 
                           soupQuantile,
                           adjustMethod,
                           roundToInt,
                           pCut) {
    # Building SoupX's SoupChannel object
    toc <- expData(inSCE, useAssay)
    if (is.null(background)) {
        scNoDrops <- SoupX::SoupChannel(toc, toc, calcSoupProfile = FALSE)
        soupProf <- data.frame(row.names = rownames(toc), 
                               est = rowSums(toc)/sum(toc), 
                               counts = rowSums(toc))
        sc <- SoupX::setSoupProfile(scNoDrops, soupProf)
    } else {
        if (inherits(background, "SingleCellExperiment")) {
            if (is.null(bgAssayName)) {
                bgAssayName <- useAssay
            }
            tod <- SummarizedExperiment::assay(background, bgAssayName)
        } else {
            tod <- background
        }
        sc <- SoupX::SoupChannel(tod, toc)
    }
    
    # Adding cluster info
    if (!is.null(cluster)) {
        if (is.character(cluster) & length(cluster) == 1) {
            cluster <- SummarizedExperiment::colData(inSCE)[[cluster]]
            names(cluster) <- colnames(inSCE)
            sc = SoupX::setClusters(sc, cluster)
        } else if (length(cluster) == ncol(inSCE)) {
            names(cluster) <- colnames(inSCE)
            sc = SoupX::setClusters(sc, cluster)
        } else {
            stop("Invalid cluster specification")
        }
    } else {
        message(paste0(date(), " ... Cluster info not given, "))
        message(paste0(date(), " ...   generating with Scran SNN method."))
        suppressMessages({
            inSCE <- scaterlogNormCounts(inSCE, "logcounts", useAssay)
            inSCE <- runNormalization(inSCE, "logcounts", "scaled", scale = TRUE)
            inSCE <- seuratFindHVG(inSCE, verbose = FALSE)
            inSCE <- getTopHVG(inSCE, "vst", altExp = "hvg")
            inSCE <- scaterPCA(inSCE, "scaled", "hvg")
            inSCE <- runScranSNN(inSCE, useReducedDim = "PCA", 
                                 clusterName = "SoupX_cluster")
        })
        
        sc = SoupX::setClusters(sc, setNames(inSCE$SoupX_cluster, 
                                             colnames(inSCE)))
    }
    
    sc <- SoupX::autoEstCont(sc, doPlot = FALSE, tfidfMin = tfidfMin, 
                             soupQuantile = soupQuantile)
    out <- SoupX::adjustCounts(sc, method = adjustMethod, 
                               roundToInt = roundToInt, pCut = pCut)
    return(list(sc = sc, out = out))
}

# TODO
# S4 method soupXResults(sce, “sample”) <- fit

# plotSoupXResults
# plotSoupXBGResults

