## Original code from DoubletFinder package
## (https://github.com/chris-mcginnis-ucsf/DoubletFinder)

.parallel_paramSweep <- function(n, n.real.cells, real.cells, pK, pN, data,
                                   orig.commands, PCs, sct, verbose, seed) {
    sweep.res.list <- list()
    list.ind <- 0

    ## Make merged real-artifical data
    if (verbose) {
        print(paste("Creating artificial doublets for pN = ",
                    pN[n] * 100, "%", sep = ""))
    }
    n_doublets <- ceiling(n.real.cells / (1 - pN[n]) - n.real.cells)
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Pre-process Seurat object
    if (sct == FALSE) {
        if (verbose) {
            print("Creating Seurat object...")
        }
        seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)

        if (verbose) {
            print("Normalizing Seurat object...")
        }
        seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets,
                                       normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                       scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                       margin = orig.commands$NormalizeData.RNA@params$margin,
                                       verbose = verbose
        )

        if (verbose) {
            print("Finding variable genes...")
        }
        seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets,
                                              selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                              loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                              clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                              mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                              dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                              num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                              binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                              nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                              mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                              dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff,
                                              verbose = verbose
        )

        if (verbose) {
            print("Scaling data...")
        }
        seu_wdoublets <- Seurat::ScaleData(
            seu_wdoublets,
            features = orig.commands$ScaleData.RNA$features,
            model.use = orig.commands$ScaleData.RNA$model.use,
            do.scale = orig.commands$ScaleData.RNA$do.scale,
            do.center = orig.commands$ScaleData.RNA$do.center,
            scale.max = orig.commands$ScaleData.RNA$scale.max,
            block.size = orig.commands$ScaleData.RNA$block.size,
            min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block,
            verbose = verbose
        )

        if (verbose) {
            print("Running PCA...")
        }
        seu_wdoublets <- Seurat::RunPCA(
            seu_wdoublets,
            features = orig.commands$ScaleData.RNA$features,
            npcs = length(PCs),
            rev.pca = orig.commands$RunPCA.RNA$rev.pca,
            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
            seed.use = seed,
            verbose = FALSE
        )
    }

    if (sct == TRUE) {
        if (verbose) {
            print("Creating Seurat object...")
        }
        seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)

        if (verbose) {
            print("Running SCTransform...")
        }
        seu_wdoublets <- Seurat::SCTransform(seu_wdoublets)

        if (verbose) {
            print("Running PCA...")
        }
        seu_wdoublets <- Seurat::RunPCA(seu_wdoublets,
                                npcs = length(PCs), seed.use = seed,
                                verbose = verbose
        )
    }

    ## Compute PC distance matrix
    if (verbose) {
        print("Calculating PC distance matrix...")
    }
    nCells <- nrow(seu_wdoublets@meta.data)
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
    seu_wdoublets <- NULL
    gc()
    dist.mat <- fields::rdist(pca.coord)[, 1:n.real.cells]

    ## Pre-order PC distance matrix prior to iterating across pK
    ## for pANN computations
    if (verbose) {
        print("Defining neighborhoods...")
    }

    for (i in 1:n.real.cells) {
        dist.mat[, i] <- order(dist.mat[, i])
    }

    ## Trim PC distance matrix for faster manipulations
    ind <- round(nCells * max(pK)) + 5
    dist.mat <- dist.mat[1:ind, ]

    ## Compute pANN across pK sweep

    if (verbose) {
        print("Computing pANN across all pK...")
    }
    for (k in 1:length(pK)) {
        if (verbose) {
            print(paste("pK = ", pK[k], "...", sep = ""))
        }
        pk.temp <- round(nCells * pK[k])
        pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
        colnames(pANN) <- "pANN"
        rownames(pANN) <- real.cells
        list.ind <- list.ind + 1

        for (i in 1:n.real.cells) {
            neighbors <- dist.mat[2:(pk.temp + 1), i]
            pANN$pANN[i] <- length(which(neighbors > n.real.cells)) / pk.temp
        }

        sweep.res.list[[list.ind]] <- pANN
    }

    return(sweep.res.list)
}

.paramSweep <- function(seu, PCs = 1:10, sct = FALSE,
                          verbose = FALSE, num.cores, seed) {
    ## Set pN-pK param sweep ranges
    pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    ## Remove pK values with too few cells
    min.cells <- round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
    pK.test <- round(pK * min.cells)
    pK <- pK[which(pK.test >= 1)]

    ## Extract pre-processing parameters from original data analysis workflow
    orig.commands <- seu@commands

    ## Down-sample cells to 10000 (when applicable) for computational effiency
    if (nrow(seu@meta.data) > 10000) {
        real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data),
                                                     10000, replace = FALSE)]
        data <- seu@assays$RNA@counts[, real.cells]
        n.real.cells <- ncol(data)
    }

    if (nrow(seu@meta.data) <= 10000) {
        real.cells <- rownames(seu@meta.data)
        data <- seu@assays$RNA@counts
        n.real.cells <- ncol(data)
    }
    ## Iterate through pN, computing pANN vectors at varying pK
    # no_cores <- detectCores()-1
    if (is.null(num.cores)) {
        num.cores <- 1
    }

    if (num.cores > 1) {
        cl <- parallel::makeCluster(num.cores)

        output2 <- withr::with_seed(
            seed,
            parallel::mclapply(as.list(1:length(pN)),
                     FUN = .parallel_paramSweep,
                     n.real.cells,
                     real.cells,
                     pK,
                     pN,
                     data,
                     orig.commands,
                     PCs,
                     sct, mc.cores = num.cores,
                     verbose = verbose,
                     seed = seed,
                     mc.set.seed = FALSE
            )
        )
        parallel::stopCluster(cl)
    } else {
        output2 <- lapply(as.list(1:length(pN)),
                          FUN = .parallel_paramSweep,
                          n.real.cells,
                          real.cells,
                          pK,
                          pN,
                          data,
                          orig.commands,
                          PCs,
                          sct,
                          seed = seed,
                          verbose = verbose
        )
    }

    ## Write parallelized output into list
    sweep.res.list <- list()
    list.ind <- 0
    for (i in 1:length(output2)) {
        for (j in 1:length(output2[[i]])) {
            list.ind <- list.ind + 1
            sweep.res.list[[list.ind]] <- output2[[i]][[j]]
        }
    }
    ## Assign names to list of results
    name.vec <- NULL
    for (j in 1:length(pN)) {
        name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
    }
    names(sweep.res.list) <- name.vec
    return(sweep.res.list)

}

.runDoubletFinder <- function(counts, seuratPcs, seuratRes, formationRate,
                              seuratNfeatures, verbose = FALSE,
                              nCores = NULL, seed = 12345) {

  ## Convert to sparse matrix if not already in that format
  counts <- .convertToMatrix(counts)
  colnames(counts) <- gsub("_", "-", colnames(counts))

  seurat <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = counts,
    project = "seurat", min.features = 0
  ))
  seurat <- Seurat::NormalizeData(
    object = seurat,
    normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = verbose
  )
  seurat <- Seurat::FindVariableFeatures(seurat,
    selection.method = "vst",
    nfeatures = seuratNfeatures,
    verbose = verbose
  )

  allGenes <- rownames(seurat)
  seurat <- Seurat::ScaleData(seurat, features = allGenes, verbose = verbose)

  numPc <- min(ncol(seurat@assays$RNA@scale.data) - 1, 50)
  seurat <- Seurat::RunPCA(seurat,
    features =
      Seurat::VariableFeatures(object = seurat),
    npcs = numPc, verbose = verbose, seed.use = seed
  )
  seurat <- Seurat::FindNeighbors(seurat, dims = seuratPcs, verbose = verbose)
  seurat <- Seurat::FindClusters(seurat,
    resolution = seuratRes,
    verbose = verbose,
    random.seed = seed
  )
  invisible(sweepResListSeurat <- .paramSweep(seurat,
    PCs = seuratPcs,
    sct = FALSE,
    num.cores = nCores,
    verbose = verbose,
    seed = seed
  ))
  sweepStatsSeurat <- .summarizeSweep(sweepResListSeurat,
    GT = FALSE
  )
  bcmvnSeurat <- .find.pK(sweepStatsSeurat)
  pkOptimal <- as.numeric(as.matrix(bcmvnSeurat$pK[
    which.max(bcmvnSeurat$MeanBC)
  ]))
  annotations <- seurat@meta.data$seurat_clusters
  homotypicProp <- .modelHomotypic(annotations)
  nExpPoi <- round(formationRate * ncol(seurat@assays$RNA))
  seurat <- .doubletFinder_v3(seurat,
    PCs = seuratPcs,
    pN = 0.25,
    pK = pkOptimal,
    nExp = nExpPoi,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  names(seurat@meta.data)[6] <- "doubletFinderAnnScore"
  names(seurat@meta.data)[7] <- "doubletFinderLabel"
  return(seurat)
}

#' @title Generates a doublet score for each cell via doubletFinder
#' @description Uses doubletFinder to determine cells within the dataset
#'  suspected to be doublets.
#' @param inSCE Input SingleCellExperiment object. Must contain a counts matrix
#' @param useAssay  Indicate which assay to use. Default "counts".
#' @param sample Numeric vector. Each cell will be assigned a sample number.
#' @param seed Seed for the random number generator. Default 12345.
#' @param seuratNfeatures Integer. Number of highly variable genes to use.
#'  Default 2000.
#' @param seuratPcs Numeric vector. The PCs used in seurat function to
#'   determine number of clusters. Default 1:15.
#' @param seuratRes Numeric vector. The resolution parameter used in seurat,
#'  which adjusts the number of clusters determined via the algorithm.
#'  Default 1.5.
#' @param formationRate Doublet formation rate used within algorithm.
#'  Default 0.075.
#' @param nCores Number of cores used for running the function.
#' @param verbose Boolean. Wheter to print messages from Seurat and
#'  DoubletFinder. Default FALSE.
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- runDoubletFinder(sce)
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment reducedDim<-
runDoubletFinder <- function(inSCE,
                             useAssay = "counts",
                             sample = NULL,
                             seed = 12345,
                             seuratNfeatures = 2000,
                             seuratPcs = 1:15,
                             seuratRes = 1.5,
                             formationRate = 0.075,
                             nCores = NULL,
                             verbose = FALSE) {

  # argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()), sys.frame(sys.nframe()))

  if (!is.null(sample)) {
    if (length(sample) != ncol(inSCE)) {
      stop("'sample' must be the same length as the number of columns in 'inSCE'")
    }
  } else {
    sample <- rep(1, ncol(inSCE))
  }

  message(paste0(date(), " ... Running 'doubletFinder'"))

  doubletScore <- rep(NA, ncol(inSCE))
  doubletLabel <- rep(NA, ncol(inSCE))
  allSampleNumbers <- sort(unique(sample))

  for (res in seuratRes) {
    output <- S4Vectors::DataFrame(
      row.names = colnames(inSCE),
      doubletFinder_doublet_score = numeric(ncol(inSCE)),
      doubletFinder_doublet_label = numeric(ncol(inSCE))
    )
    umapDims <- matrix(ncol = 2,
                       nrow = ncol(inSCE))
    rownames(umapDims) = colnames(inSCE)
    colnames(umapDims) = c("UMAP_1", "UMAP_2")

    ## Loop through each sample and run doubletFinder
    samples <- unique(sample)

    for (i in seq_len(length(samples))) {
      sceSampleInd <- sample == samples[i]
      sceSample <- inSCE[, sceSampleInd]
      mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
      if(length(seuratPcs) > ncol(mat)){
        seuratPcs <- 1:ncol(mat)
      }

      result <- suppressMessages(withr::with_seed(
        seed,
        .runDoubletFinder(
          counts = mat,
          seuratPcs = seuratPcs,
          seuratRes = res,
          seuratNfeatures = seuratNfeatures,
          formationRate = formationRate,
          nCores = nCores,
          verbose = verbose,
          seed = seed
        )
      ))
      result <- suppressMessages(Seurat::RunUMAP(result,
                                                 dims = 1:10,
                                                 verbose = verbose,
                                                 umap.method = "umap-learn",
                                                 metric = "correlation"))

      seuratDims <- Seurat::Embeddings(result, reduction = "umap")
      umapDims[sceSampleInd, 1] <- seuratDims[,1]
      umapDims[sceSampleInd, 2] <- seuratDims[,2]
      output[sceSampleInd, 1] <- result@meta.data$doubletFinderAnnScore
      output[sceSampleInd, 2] <- result@meta.data$doubletFinderLabel
    }

    colnames(output) <- paste0(colnames(output), "_resolution_", res)

    argsList <- argsList[!names(argsList) %in% ("...")]
    inSCE@metadata$runDoubletFinder <- argsList[-1]
    inSCE@metadata$runDoubletFinder$packageVersion <- "2.0.2"
    colData(inSCE) <- cbind(colData(inSCE), output)
    reducedDim(inSCE,'doubletFinder_UMAP') <- umapDims
  }
  return(inSCE)
}


.summarizeSweep <- function(sweep.list, GT = FALSE, GT.calls = NULL) {
  #require(KernSmooth); require(ROCR)
  ## Set pN-pK param sweep ranges
  name.vec <- names(sweep.list)
  name.vec <- unlist(strsplit(name.vec, split="pN_"))
  name.vec <- name.vec[seq(2, length(name.vec), by=2)]
  name.vec <- unlist(strsplit(name.vec, split="_pK_"))
  pN <- as.numeric(unique(name.vec[seq(1, length(name.vec), by=2)]))
  pK <- as.numeric(unique(name.vec[seq(2, length(name.vec), by=2)]))

  ## Initialize data structure w/ or w/o AUC column, depending on whether ground-truth doublet classifications are available
  if (GT == TRUE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=4))
    colnames(sweep.stats) <- c("pN","pK","AUC","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }

  if (GT == FALSE) {
    sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=3))
    colnames(sweep.stats) <- c("pN","pK","BCreal")
    sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
    sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
  }

  ## Perform pN-pK parameter sweep summary
  for (i in 1:length(sweep.list)) {
    res.temp <- sweep.list[[i]]

    ## Use gaussian kernel density estimation of pANN vector to compute bimodality coefficient
    gkde <- stats::approxfun(KernSmooth::bkde(res.temp$pANN, kernel="normal"))
    x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
    sweep.stats$BCreal[i] <- .bimodality_coefficient(gkde(x))

    if (GT == FALSE) { next }

    ## If ground-truth doublet classifications are available, perform ROC analysis on logistic
    ## regression model trained using pANN vector
    meta <- as.data.frame(matrix(0L, nrow=nrow(res.temp), ncol=2))
    meta[,1] <- GT.calls
    meta[,2] <- res.temp$pANN
    train.ind <- sample(1:nrow(meta), round(nrow(meta)/2), replace=FALSE)
    test.ind <- (1:nrow(meta))[-train.ind]
    colnames(meta) <- c("SinDub","pANN")
    meta$SinDub <- factor(meta$SinDub, levels = c("Doublet","Singlet"))
    model.lm <- stats::glm(SinDub ~ pANN, family=stats::binomial(link='logit'), data=meta, subset=train.ind)
    prob <- stats::predict(model.lm, newdata=meta[test.ind, ], type="response")
    ROCpred <- ROCR::prediction(predictions=prob, labels=meta$SinDub[test.ind])
    perf.auc <- ROCR::performance(ROCpred, measure="auc")
    sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
  }

  return(sweep.stats)
}


.bimodality_coefficient <- function(x) {
  G <- .skewness(x)
  sample.excess.kurtosis <- .kurtosis(x)
  K <- sample.excess.kurtosis
  n <- length(x)
  B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
  return(B)
}

.skewness <- function(x) {
  n <- length(x)
  S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  S <- S*(sqrt(n*(n-1)))/(n-2)
  return(S)
}

.kurtosis <- function (x) {
  n <- length(x)
  K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
  K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
  return(K)
}

.modelHomotypic <- function(annotations) {
  anno.freq <- table(annotations)/length(annotations)
  x <- sum(anno.freq^2)
  return(x)
}

.find.pK <- function(sweep.stats) {

  ## Implementation for data without ground-truth doublet classifications
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- stats::sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(stats::sd(sweep.stats[ind, "BCreal"])^2)
    }
    return(bc.mvn)
  }

  ## Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
    colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- stats::sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(stats::sd(sweep.stats[ind, "BCreal"])^2)
    }

    return(bc.mvn)

  }
}

.doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE) {

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Store important pre-processing information
    orig.commands <- seu@commands

    ## Pre-process Seurat object
    if (sct == FALSE) {
      seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)

      seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)

      seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

      seu_wdoublets <- Seurat::ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

      seu_wdoublets <- Seurat::RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }

    if (sct == TRUE) {
      #require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- Seurat::SCTransform(seu_wdoublets)

      print("Running PCA...")
      seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }

    ## Compute PC distance matrix
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
}
