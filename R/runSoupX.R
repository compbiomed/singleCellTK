#' @title Detecting and correct contamination with SoupX
#' @description A wrapper function for \link[SoupX]{autoEstCont} and 
#' \link[SoupX]{adjustCounts}. Identify potential contamination from 
#' experimental factors such as ambient RNA. Visit
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
#' When not supplied, \code{\link[scran]{quickCluster}} method will be applied.
#' @param reducedDimName A single character string of the prefix of output 
#' corrected embedding matrix for each sample. Default \code{"SoupX_UMAP_"} when 
#' not using a background, otherwise, \code{"SoupX_bg_UMAP_"}.
#' @param tfidfMin Numeric. Minimum value of tfidf to accept for a marker gene. 
#' Default \code{1}. See \code{?SoupX::autoEstCont}. 
#' @param soupQuantile Numeric. Only use genes that are at or above this 
#' expression quantile in the soup. This prevents inaccurate estimates due to 
#' using genes with poorly constrained contribution to the background. Default 
#' \code{0.9}. See \code{?SoupX::autoEstCont}. 
#' @param maxMarkers Integer. If we have heaps of good markers, keep only the 
#' best maxMarkers of them. Default \code{100}. See \code{?SoupX::autoEstCont}. 
#' @param contaminationRange Numeric vector of two elements. This constrains 
#' the contamination fraction to lie within this range. Must be between 0 and 1.
#' The high end of this range is passed to 
#' \code{\link[SoupX]{estimateNonExpressingCells}} as 
#' \code{maximumContamination}. Default \code{c(0.01, 0.8)}. See 
#' \code{?SoupX::autoEstCont}. 
#' @param rhoMaxFDR Numeric. False discovery rate passed to 
#' \code{\link[SoupX]{estimateNonExpressingCells}}, to test if rho is less than 
#' \code{maximumContamination}. Default \code{0.2}. See 
#' \code{?SoupX::autoEstCont}. 
#' @param priorRho Numeric. Mode of gamma distribution prior on contamination 
#' fraction. Default \code{0.05}. See \code{?SoupX::autoEstCont}. 
#' @param priorRhoStdDev Numeric. Standard deviation of gamma distribution prior 
#' on contamination fraction. Default \code{0.1}. See 
#' \code{?SoupX::autoEstCont}. 
#' @param forceAccept Logical. Should we allow very high contamination fractions
#' to be used. Passed to \code{\link[SoupX]{setContaminationFraction}}. Default 
#' \code{FALSE}. See \code{?SoupX::autoEstCont}. 
#' @param adjustMethod Character. Method to use for correction. One of 
#' \code{'subtraction'}, \code{'soupOnly'}, or \code{'multinomial'}. Default 
#' \code{'subtraction'}. See \code{?SoupX::adjustCounts}. 
#' @param roundToInt Logical. Should the resulting matrix be rounded to 
#' integers? Default \code{FALSE}. See \code{?SoupX::adjustCounts}. 
#' @param tol Numeric. Allowed deviation from expected number of soup counts. 
#' Don't change this. Default \code{0.001}. See \code{?SoupX::adjustCounts}. 
#' @param pCut Numeric. The p-value cut-off used when 
#' \code{method = 'soupOnly'}. Default \code{0.01}. See 
#' \code{?SoupX::adjustCounts}. 
#' @return The input \code{inSCE} object with \code{soupX_nUMIs}, 
#' \code{soupX_clustrers}, \code{soupX_contamination} appended to \code{colData} 
#' slot; \code{soupX_{sample}_est} and \code{soupX_{sample}_counts} for each 
#' sample appended to \code{rowData} slot; and other computational metrics at 
#' \code{getSoupX(inSCE)}. Replace "soupX" to "soupX_bg" when \code{background}
#' is used. 
#' @export
#' @author Yichen Wang
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' # SoupX does not work for toy example, 
#' # can be tested with `sce <- importExampleData("pbmc3k")`
#' sce <- runSoupX(sce, sample = "sample")
#' }
runSoupX <- function(inSCE, 
                     sample = NULL,
                     useAssay = "counts", 
                     background = NULL, 
                     bgAssayName = NULL,
                     bgBatch = NULL,
                     assayName = ifelse(is.null(background), 
                                        "SoupX", "SoupX_bg"), 
                     cluster = NULL, 
                     reducedDimName = ifelse(is.null(background), 
                                        "SoupX_UMAP_", "SoupX_bg_UMAP_"), 
                     tfidfMin = 1, 
                     soupQuantile = 0.9,
                     maxMarkers = 100,
                     contaminationRange = c(0.01, 0.8),
                     rhoMaxFDR = 0.2,
                     priorRho = 0.05,
                     priorRhoStdDev = 0.1,
                     forceAccept = FALSE,
                     adjustMethod = c("subtraction", "soupOnly", "multinomial"),
                     roundToInt = FALSE,
                     tol = 0.001,
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
                if (!all(sample %in% bgBatch)) {
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
                                  reducedDimName = reducedDimName,
                                  tfidfMin = tfidfMin,
                                  soupQuantile = soupQuantile,
                                  maxMarkers = maxMarkers,
                                  contaminationRange = contaminationRange,
                                  rhoMaxFDR = rhoMaxFDR,
                                  priorRho = priorRho,
                                  priorRhoStdDev = priorRhoStdDev,
                                  forceAccept = forceAccept,
                                  adjustMethod = adjustMethod,
                                  roundToInt = roundToInt,
                                  tol = tol,
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
                                  reducedDimName = reducedDimName,
                                  tfidfMin = tfidfMin,
                                  soupQuantile = soupQuantile,
                                  maxMarkers = maxMarkers,
                                  contaminationRange = contaminationRange,
                                  rhoMaxFDR = rhoMaxFDR,
                                  priorRho = priorRho,
                                  priorRhoStdDev = priorRhoStdDev,
                                  forceAccept = forceAccept,
                                  adjustMethod = adjustMethod,
                                  roundToInt = roundToInt,
                                  tol = tol,
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
        meta <- results[[s]]$sc$fit
        meta$nDropUMIs <- results[[s]]$sc$nDropUMIs
        meta$param <- list(sample = sample, useAssay = useAssay, 
                           bgAssayName = bgAssayName, bgBatch = bgBatch,
                           assayName = assayName, 
                           tfidfMin = tfidfMin,
                           soupQuantile = soupQuantile,
                           maxMarkers = maxMarkers,
                           contaminationRange = contaminationRange,
                           rhoMaxFDR = rhoMaxFDR,
                           priorRho = priorRho,
                           priorRhoStdDev = priorRhoStdDev,
                           forceAccept = forceAccept,
                           adjustMethod = adjustMethod,
                           roundToInt = roundToInt,
                           tol = tol,
                           pCut = pCut,
                           reducedDimName = paste0(reducedDimName, s),
                           sessionInfo = utils::sessionInfo())
        # Output inSCE need to have separated UMAP calculated for each sample
        sampleUMAP <- matrix(nrow = ncol(inSCE), ncol = 2)
        rownames(sampleUMAP) <- colnames(inSCE)
        sampleUMAP[sampleIdx[[s]]] <- results[[s]]$umap
        if (!is.null(cluster)) {
            meta$param$cluster <- cluster
        } else {
            if (!is.null(background)) {
                meta$param$cluster <- "soupX_bg_clusters"
            } else {
                meta$param$cluster <- "soupX_clusters"
            }
        }
        if (!is.null(background)) {
            newColData$soupX_bg_nUMIs[sampleIdx[[s]]] <- results[[s]]$sc$metaData$nUMIs
            newColData$soupX_bg_clusters[sampleIdx[[s]]] <- paste0(s, "-", results[[s]]$sc$metaData$clusters)
            newColData$soupX_bg_contamination[sampleIdx[[s]]] <- results[[s]]$sc$metaData$rho
            newRowData[[paste0("soupX_bg_",s,"_est")]] <- results[[s]]$sc$soupProfile$est
            newRowData[[paste0("soupX_bg_",s,"_counts")]] <- results[[s]]$sc$soupProfile$counts
            getSoupX(inSCE, sampleID = s, background = TRUE) <- meta
        } else {
            newColData$soupX_nUMIs[sampleIdx[[s]]] <- results[[s]]$sc$metaData$nUMIs
            newColData$soupX_clusters[sampleIdx[[s]]] <- paste0(s, "-", results[[s]]$sc$metaData$clusters)
            newColData$soupX_contamination[sampleIdx[[s]]] <- results[[s]]$sc$metaData$rho
            newRowData[[paste0("soupX_",s,"_est")]] <- results[[s]]$sc$soupProfile$est
            newRowData[[paste0("soupX_",s,"_counts")]] <- results[[s]]$sc$soupProfile$counts
            getSoupX(inSCE, sampleID = s) <- meta
        }
        SingleCellExperiment::reducedDim(inSCE, 
                                         paste0(reducedDimName, 
                                                s)) <- sampleUMAP
    }
    expData(inSCE, assayName, tag = "raw") <- corrAssay
    inSCE <- expSetDataTag(inSCE, "raw", assayName)
    SummarizedExperiment::colData(inSCE) <- newColData
    SummarizedExperiment::rowData(inSCE) <- newRowData
    return(inSCE)
}

.SoupXOneBatch <- function(inSCE, 
                           useAssay, 
                           background, 
                           bgAssayName,
                           cluster,
                           reducedDimName,
                           tfidfMin, 
                           soupQuantile,
                           maxMarkers,
                           contaminationRange,
                           rhoMaxFDR,
                           priorRho,
                           priorRhoStdDev,
                           forceAccept,
                           adjustMethod,
                           roundToInt,
                           tol,
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
            sc <- SoupX::setClusters(sc, cluster)
        } else if (length(cluster) == ncol(inSCE)) {
            names(cluster) <- colnames(inSCE)
            sc <- SoupX::setClusters(sc, cluster)
        } else {
            stop("Invalid cluster specification")
        }
    } else {
        message(paste0(date(), " ... Cluster info not given, "))
        message(paste0(date(), " ...   generating with Scran SNN method."))
        suppressMessages({
            c <- scran::quickCluster(inSCE, assay.type = useAssay, 
                                     method = "igraph")
            inSCE$SoupX_cluster <- c
        })
        sc <- SoupX::setClusters(sc, stats::setNames(inSCE$SoupX_cluster, 
                                                    colnames(inSCE)))
    }
    
    sc <- SoupX::autoEstCont(sc, doPlot = FALSE, tfidfMin = tfidfMin, 
                             soupQuantile = soupQuantile, 
                             maxMarkers = maxMarkers,
                             contaminationRange = contaminationRange,
                             rhoMaxFDR = rhoMaxFDR, priorRho = priorRho,
                             priorRhoStdDev = priorRhoStdDev,
                             forceAccept = forceAccept)
    
    out <- SoupX::adjustCounts(sc, method = adjustMethod, 
                               roundToInt = roundToInt, tol = tol, pCut = pCut)
    message(paste0(date(), " ... Generating UMAP"))
    inSCE <- getUMAP(inSCE, useAssay = useAssay, reducedDimName = "sampleUMAP")
    return(list(sc = sc, out = out, 
                umap = SingleCellExperiment::reducedDim(inSCE, "sampleUMAP")))
}

#' @title Get or Set SoupX Result
#' @rdname getSoupX
#' @description S4 method for getting and setting SoupX results that cannot be
#' appended to either \code{rowData(inSCE)} or \code{colData(inSCE)}. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object. For getter method,
#' \code{\link{runSoupX}} must have been already applied. 
#' @param sampleID Character vector. For getter method, the samples that should 
#' be included in the returned list. Leave this \code{NULL} for all samples. 
#' Default \code{NULL}. For setter method, only one sample allowed. 
#' @param background Logical. Whether \code{background} was applied when 
#' running \code{\link{runSoupX}}. Default \code{FALSE}.
#' @param value Dedicated list object of SoupX results. 
#' @return For getter method, a list with SoupX results for specified samples. 
#' For setter method, \code{inSCE} with SoupX results updated.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' # SoupX does not work for toy example, 
#' # can be tested with `sce <- importExampleData("pbmc3k")`
#' sce <- runSoupX(sce, sample = "sample")
#' soupXResults <- getSoupX(sce)
#' }
setGeneric("getSoupX<-", function(inSCE, sampleID, background = FALSE, value) 
    standardGeneric("getSoupX<-") )

#' @title Insert SoupX result to SCE object
#' @rdname getSoupX
#' @export
setGeneric("getSoupX", function(inSCE, sampleID = NULL, background = FALSE) 
    standardGeneric("getSoupX") )

#' @title Get or Set SoupX Result
#' @rdname getSoupX
#' @description S4 method for getting and setting SoupX results that cannot be
#' appended to either \code{rowData(inSCE)} or \code{colData(inSCE)}. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object. For getter method,
#' \code{\link{runSoupX}} must have been already applied. 
#' @param sampleID Character vector. For getter method, the samples that should 
#' be included in the returned list. Leave this \code{NULL} for all samples. 
#' Default \code{NULL}. For setter method, only one sample allowed. 
#' @param background Logical. Whether \code{background} was applied when 
#' running \code{\link{runSoupX}}. Default \code{FALSE}.
#' @param value Dedicated list object of SoupX results. 
#' @return For getter method, a list with SoupX results for specified samples. 
#' For setter method, \code{inSCE} with SoupX results updated.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' # SoupX does not work for toy example, 
#' # can be tested with `sce <- importExampleData("pbmc3k")`
#' sce <- runSoupX(sce, sample = "sample")
#' soupXResults <- getSoupX(sce)
#' }
setMethod("getSoupX", 
          "SingleCellExperiment", 
          function(inSCE, 
                   sampleID = NULL,
                   background = FALSE){
    if (isTRUE(background)) {
        all.results <- S4Vectors::metadata(inSCE)$sctk$SoupX_bg
    } else {
        all.results <- S4Vectors::metadata(inSCE)$sctk$SoupX
    }
    if(is.null(all.results)) {
        stop("No result from 'SoupX' is found. Please run `runSoupX()` first, ",
             "or check the setting of `background`.") 
    }
    results <- all.results
    if (!is.null(sampleID)) {
        if (!all(sampleID  %in% names(all.results))) {
            stop("Sample(s) not found in the results for tool 'SoupX': ",
                 paste(sampleID[!sampleID %in% names(all.results)], collapse = ", "))
        }
        results <- all.results[sampleID]
    }
    return(results)
})

#' @title Insert SoupX result to SCE object
#' @rdname getSoupX
#' @export
setReplaceMethod("getSoupX", 
                 c("SingleCellExperiment"), 
                 function(inSCE, 
                          sampleID, 
                          background = FALSE, 
                          value) {
                     if (isTRUE(background)) {
                         inSCE@metadata$sctk$SoupX_bg[[sampleID]] <- value
                     } else {
                         inSCE@metadata$sctk$SoupX[[sampleID]] <- value
                     }
                     return(inSCE)
                 })

#' Plot SoupX Result
#' @description This function will generate a combination of plots basing on the 
#' correction done by SoupX. For each sample, there will be a UMAP with cluster 
#' labeling, followed by a number of UMAPs showing the change in selected top 
#' markers. The cluster labeling is what should be used for SoupX to estimate 
#' the contamination. The Soup Fraction is calculated by subtracting the gene 
#' expression value of the output corrected matrix from that of the original 
#' input matrix, and then devided by the input. 
#' @param inSCE A \linkS4class{SingleCellExperiment} object. With 
#' \code{\link{runSoupX}} already applied. 
#' @param sample Character vector. Indicates which sample each cell belongs to. 
#' Default \code{NULL}.
#' @param background Logical. Whether \code{background} was applied when 
#' running \code{\link{runSoupX}}. Default \code{FALSE}.
#' @param reducedDimName Character. The embedding to use for plotting. Leave it
#' \code{NULL} for using the sample-specific UMAPs generated when running 
#' \code{\link{runSoupX}}. Default \code{NULL}. 
#' @param plotNCols Integer. Number of columns for the plot grid per sample. 
#' Will determine the number of top markers to show together with 
#' \code{plotNRows}. Default \code{3}.
#' @param plotNRows Integer. Number of rows for the plot grid per sample. Will
#' determine the number of top markers to show together with \code{plotNCols}.
#' Default \code{2}.
#' @param baseSize Numeric. The base font size for all text. Default 12. Can be 
#' overwritten by titleSize, axisSize, and axisLabelSize, legendSize, 
#' legendTitleSize. Default \code{8}.
#' @param combinePlot Must be either \code{"all"}, \code{"sample"}, or 
#' \code{"none"}. \code{"all"} will combine all plots into a single 
#' \code{.ggplot} object, while \code{"sample"} will output a list of plots 
#' separated by sample. Default \code{"all"}.
#' @param xlab Character vector. Label for x-axis. Default \code{NULL}.
#' @param ylab Character vector. Label for y-axis. Default \code{NULL}.
#' @param dim1 See \code{\link{plotSCEDimReduceColData}}. Default \code{NULL}.
#' @param dim2 See \code{\link{plotSCEDimReduceColData}}. Default \code{NULL}.
#' @param labelClusters Logical. Whether the cluster labels are plotted. Default
#' \code{FALSE}.
#' @param clusterLabelSize Numeric. Determines the size of cluster label when 
#' \code{labelClusters} is set to \code{TRUE}. Default \code{3.5}.
#' @param defaultTheme Logical. Adds grid to plot when \code{TRUE}. Default 
#' \code{TRUE}.
#' @param dotSize Numeric. Size of dots. Default \code{0.5}.
#' @param transparency Numeric. Transparency of the dots, values will be from 0
#' to 1. Default \code{1}.
#' @param titleSize Numeric. Size of title of plot. Default \code{15}.
#' @param axisLabelSize Numeric. Size of x/y-axis labels. Default \code{NULL}.
#' @param axisSize Numeric. Size of x/y-axis ticks. Default \code{NULL}.
#' @param legendSize Numeric. Size of legend. Default \code{NULL}.
#' @param legendTitleSize Numeric. Size of legend title. Default \code{NULL}.
#' @return ggplot object of the combination of UMAPs. See description.
#' @export
#' @examples 
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' # SoupX does not work for toy example, 
#' # can be tested with `sce <- importExampleData("pbmc3k")`
#' sce <- runSoupX(sce, sample = "sample")
#' plotSoupXResults(sce)
#' }
plotSoupXResults <- function(inSCE, 
                            sample = NULL, 
                            background = FALSE, 
                            reducedDimName=NULL,
                            plotNCols = 3,
                            plotNRows = 2,
                            baseSize=8,
                            combinePlot = c("all", "sample", "none"),
                            xlab=NULL,
                            ylab=NULL,
                            dim1=NULL,
                            dim2=NULL,
                            labelClusters = FALSE,
                            clusterLabelSize = 3.5,
                            defaultTheme=TRUE,
                            dotSize=0.5,
                            transparency=1,
                            titleSize=NULL,
                            axisLabelSize=NULL,
                            axisSize=NULL,
                            legendSize=NULL,
                            legendTitleSize=NULL
                            )
{
    combinePlot <- match.arg(combinePlot)
    sampleID <- unique(sample)
    results <- getSoupX(inSCE, sampleID = sampleID, background = background)
    # Doing this redundancy-like step because: If sample given NULL when running
    # runSoupX(), the actual sample label saved will be "all_cells", which users
    # won't know. 
    samples <- names(results)
    samplePlots <- list()
    for (s in samples) {
        sampleRes <- results[[s]]
        param <- sampleRes$param
        sampleIdx <- param$sample == s
        tmpSCE <- inSCE[,sampleIdx]
        markerTable <- results[[s]]$markersUsed
        markerTable <- markerTable[order(markerTable$qval),]
        # Iteratively pop out the top significant marker from each cluster.
        # i.e. Get top marker from c1 to cn, then the second top marker again
        # from c1 to cn, as long as there is still a marker for cn. 
        uniqCluster <- unique(markerTable$cluster)
        markerToUse <- character()
        markerClusterMap <- character()
        nMarkerToUse <- plotNCols * plotNRows - 1
        nIter <- 0
        while (length(markerToUse) < nMarkerToUse) {
            nIter <- nIter + 1
            # c is the cluster to look at in this iteration
            c <- uniqCluster[(nIter-1)%%length(uniqCluster) + 1]
            markerPerCluster <- markerTable[markerTable$cluster == c,]
            topMarker <- markerPerCluster[1, "gene"]
            markerToUse <- c(markerToUse, topMarker)
            markerClusterMap <- c(markerClusterMap, c)
            markerTable <- markerTable[-which(markerTable$gene == topMarker),]
        }
        names(markerClusterMap) <- markerToUse
        # Get the per-sample UMAP if not specifying one
        useRedDim <- reducedDimName
        if (is.null(reducedDimName)) {
            useRedDim <- param$reducedDimName
        }
        plotList <- list(
            scatter_soupXClusters = 
                plotSCEDimReduceColData(tmpSCE, 
                                        param$cluster, 
                                        useRedDim,
                                        title = "Cluster",
                                        labelClusters = labelClusters,
                                        clusterLabelSize = clusterLabelSize,
                                        xlab=xlab, 
                                        ylab=ylab,
                                        dim1=dim1,
                                        dim2=dim2,
                                        defaultTheme=defaultTheme,
                                        dotSize=dotSize,
                                        transparency=transparency,
                                        baseSize=baseSize,
                                        titleSize=titleSize,
                                        axisLabelSize=axisLabelSize,
                                        axisSize=axisSize,
                                        legendSize=legendSize,
                                        legendTitleSize=legendTitleSize)
            )
        for (g in markerToUse) {
            # Soup fraction was calculated basing on SoupX's original method.
            # Credit to SoupX::plotChangeMap
            oldName <- param$useAssay
            newName <- param$assayName
            old <- colSums(assay(tmpSCE, oldName)[g, , drop = FALSE])
            new <- colSums(assay(tmpSCE, newName)[g, , drop = FALSE])
            relChange = (old - new)/old
            zLims = c(0, 1)
            relChange[which(relChange < zLims[1])] = zLims[1]
            relChange[which(relChange > zLims[2])] = zLims[2]
            relChange[which(is.na(relChange))] = 0
            # Start to generate the plot
            tmpSCE$Soup_Frac <- relChange
            legendTitle <- paste0(markerClusterMap[g], ":", g, ", ", s)
            plotList[[g]] <- plotSCEDimReduceColData(tmpSCE, 
                                                     "Soup_Frac", 
                                                     useRedDim,
                                                     legendTitle = "Soup_Frac",
                                                     title = legendTitle,
                                                     xlab=xlab, 
                                                     ylab=ylab,
                                                     dim1=dim1,
                                                     dim2=dim2,
                                                     defaultTheme=defaultTheme,
                                                     dotSize=dotSize,
                                                     transparency=transparency,
                                                     baseSize=baseSize,
                                                     titleSize=titleSize,
                                                     axisLabelSize=axisLabelSize,
                                                     axisSize=axisSize,
                                                     legendSize=legendSize,
                                                     legendTitleSize=legendTitleSize
                                                     )
        }
        if (combinePlot %in% c("sample", "all")) {
            samplePlots[[s]] <- .ggSCTKCombinePlots(plotList,
                                                    plotNCols)
        } else if (combinePlot == "none") {
            samplePlots[[s]] <- plotList
        }
    }
    finalPlotList <- list(Sample = samplePlots)
    if (combinePlot == "all") {
        finalPlotList <- .ggSCTKCombinePlots(finalPlotList$Sample, ncols = 1)
    }
    if (length(samples) == 1) {
        if (combinePlot %in% c("none", "sample")) {
            finalPlotList <- finalPlotList$Sample[[1]]
        } 
    }
    return(finalPlotList)
}
