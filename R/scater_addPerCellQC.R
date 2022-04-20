#' @title Wrapper for calculating QC metrics with scater.
#' @description A wrapper function for \link[scater]{addPerCellQC}. Calculate
#'  general quality control metrics for each cell in the count matrix.
#' @param inSCE Input \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param useAssay A string specifying which assay in the SCE to use. Default
#' \code{"counts"}.
#' @param collectionName Character. Name of a \code{GeneSetCollection} obtained
#' by using one of the importGeneSet* functions. Default \code{NULL}.
#' @param geneSetList List of gene sets to be quantified. The genes in the
#' assays will be matched to the genes in the list based on
#' \code{geneSetListLocation}. Default \code{NULL}.
#' @param geneSetListLocation Character or numeric vector. If set to 'rownames',
#' then the genes in 'geneSetList' will be looked up in \code{rownames(inSCE)}.
#' If another character is supplied, then genes will be looked up in the column
#'  names of \code{rowData(inSCE)}. A character vector with the same length as
#'  \code{geneSetList} can be supplied if the IDs for different
#' gene sets are found in different places, including a mixture of 'rownames'
#' and \code{rowData(inSCE)}. An integer or integer vector can be supplied to
#' denote the column index in \code{rowData(inSCE)}. Default 'rownames'.
#' @param geneSetCollection Class of \code{GeneSetCollection} from package
#' \code{GSEAbase}. The location of the gene IDs in \code{inSCE} should be in
#' the \code{description} slot of each gene set and should follow the
#' same notation as \code{geneSetListLocation}. The function \link{getGmt} can
#' be used to read in gene sets from a GMT file. If reading a GMT file, the
#' second column for each gene set should be the description denoting the
#' location of the gene IDs in \code{inSCE}. These gene sets will be included
#' with those from \code{geneSetList} if both parameters are provided.
#' @param mitoRef Character. The species used to extract mitochondrial genes ID from 
#' build-in mitochondrial geneset in SCTK. Available species options are "human"  
#' and "mouse". Default is \code{NULL}. 
#' @param mitoIDType Character. Types of mitochondrial gene id. Now it supports "symbol", 
#' "entrez", "ensembl" and "ensemblTranscriptID". It is used with \code{mitoRef} to extract
#' mitochondrial genes from build-in mitochondrial geneset in SCTK. Default \code{NULL}. 
#' @param mitoGeneLocation Character. Describes the location within \code{inSCE} where
#' the gene identifiers in the mitochondrial gene sets should be mapped.
#' If set to \code{"rownames"} then the features will
#' be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. See \link{featureIndex} for more information.
#' Default \code{NULL}.
#' @param mitoPrefix Character. The prefix used to get mitochondrial gene from 
#' either rownames(inSCE) or columns of rowData(inSCE) specified by mitoGeneLocation. 
#' This parameter is usually used to extract mito genes from gene symbol. For example,
#' mitoPrefix = "^MT-" can be used to detect mito gene symbols like "MT-ND4".
#' @param mitoID Character. A vector of mitochondrial genes to be quantified.  
#' @param percent_top An integer vector. Each element is treated as a
#' number of top genes to compute the percentage of library size occupied by
#' the most highly expressed genes in each cell.
#' @param use_altexps Logical scalar indicating whether QC statistics should
#' be computed for alternative Experiments in x. If TRUE, statistics are
#' computed
#' for all alternative experiments.
#' Alternatively, an integer or character vector specifying the alternative
#' Experiments to use to compute QC statistics.
#' Alternatively NULL, in which case alternative experiments are not used.
#' @param flatten Logical scalar indicating whether the nested
#' \link[S4Vectors]{DataFrame-class}
#' in the output should be flattened.
#' @param detectionLimit A numeric scalar specifying the lower detection limit
#' for expression.
#' @param BPPARAM A \link{BiocParallelParam} object specifying
#' whether the QC calculations should be parallelized.
#' @details 
#' This function allows multiple ways to import mitochondrial genes and quantify 
#' their expression. 
#' \itemize{
#'   \item Using \code{mitoRef}, \code{mitoIDType} and \code{mitoGeneLocation} 
#'   parameters will load the build-in mitochondrial geneset in SCTK package. 
#'   \item Using \code{mitoPrefix} and \code{mitoGeneLocation} parameters will
#'   extract mitochondrial genes from either rownames(inSCE) or columns of 
#'   rowData(inSCE) specified ny parameter \code{mitoGeneLocation}
#'   \item Using \code{mitoID} and \code{mitoGeneLocation} parameters will quantify
#'   the expression of mitochondrial genes stored in \code{mitoID}. 
#' }
#' \code{mitoGeneLocation} is required if you use any methods mentioned above to 
#' quantify mitochondrial gene expression. Please make sure \code{mitoGeneLocation}
#' is pointing to the location within inSCE object that stores the correct mitochondrial
#' genes ID.  
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  cell QC metrics added to the \link{colData} slot. If \code{geneSetList} or
#'  \code{geneSetCollection} are provided, then the rownames for each gene set
#'  will be saved in \code{metadata(inSCE)$scater$addPerCellQC$geneSets}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' mito.ix = grep("^MT-", rowData(sce)$feature_name)
#' geneSet <- list("Mito"=rownames(sce)[mito.ix])
#' sce <- runPerCellQC(sce, geneSetList = geneSet)
#' @export
#' @importFrom SummarizedExperiment rowData colData
runPerCellQC <- function(inSCE,
                         useAssay = "counts",
                         collectionName = NULL,
                         geneSetList = NULL,
                         geneSetListLocation = "rownames",
                         geneSetCollection = NULL,
                         mitoRef = NULL,
                         mitoIDType = NULL,
                         mitoPrefix = NULL,
                         mitoID = NULL,
                         mitoGeneLocation = NULL,
                         percent_top = c(50, 100, 200, 500),
                         use_altexps = FALSE,
                         flatten = TRUE,
                         detectionLimit = 0,
                         BPPARAM = BiocParallel::SerialParam()
) {

  message(paste0(date(), " ... Running 'perCellQCMetrics'"))
  #argsList <- as.list(formals(fun = sys.function(sys.parent()),
  #                            envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Add mito gene collection from built-in mito gene sets in SCTK package
  if (is.null(mitoGeneLocation)) {
    #message("'mitoGeneLocation' not specified or specified as NULL.")
  } else {
    mitoGS <- NULL
    if (!is.null(mitoRef) & !is.null(mitoIDType)) {
      message(paste0(date(), " ... Importing Mitochondiral genes using parameter 'mitoRef', 'mitoIDType' and 'mitoGeneLocation'"))
      inSCE <- importMitoGeneSet(inSCE, reference = mitoRef, id = mitoIDType, 
                                 by = mitoGeneLocation, collectionName = "mito")
      mitoGS <- S4Vectors::metadata(inSCE)$sctk$genesets$mito[[1]]
    }

    if (mitoGeneLocation == "rownames") {
      features <- rownames(inSCE)
    } else if (mitoGeneLocation %in% names(SummarizedExperiment::rowData(inSCE))) {
      features <- SummarizedExperiment::rowData(inSCE)[[mitoGeneLocation]]
    } else {
      warning(paste0(mitoGeneLocation, " is not 'rownames' of inSCE or is not found in column names of 'rowData(inSCE)'. 
        Ignore mitoPrefix, mitoID and mitoGeneLocation parameters."))    
    }    

    mitoG <- NULL
    ## Add mito genes by greping mitoPrefix from the mitoGeneLocation
    if (!is.null(mitoPrefix)) {
      message(paste0(date(), " ... Importing Mitochondiral genes using parameter 'mitoPrefix' and 'mitoGeneLocation'"))
      mitoG <- grep(mitoPrefix, features, value = TRUE)
    }

    ## Add mito genes specified in mitoID
    if (!is.null(mitoID)) {
      message(paste0(date(), " ... Importing Mitochondiral genes using parameter 'mitoID' and 'mitoGeneLocation'"))
      mitoG <- intersect(mitoID, features)
    }

    if (!is.null(mitoG)) {
      mitoGL <- list('mito' = mitoG)
      inSCE <- importGeneSetsFromList(inSCE, geneSetList = mitoGL, 
                                      collectionName = "mito", 
                                      by = mitoGeneLocation)
      mitoGS <- S4Vectors::metadata(inSCE)$sctk$genesets$mito[[1]]
    } 

    if (is.null(mitoGS)) {
      message("No mitochondrial gene is found in 'rownames(inSCE)' or in the column of 'rowData(inSCE)' specified by ", 
        "'mitoGeneLocation'. Skip quantifying mitochondrial genes.")      
    } else {
      if (!is.null(geneSetCollection)) {
        geneSetCollection <- GSEABase::GeneSetCollection(c(mitoGS, unlist(geneSetCollection)))
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
    geneSetCollectionLocation <- vapply(geneSetCollection, GSEABase::description,
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
          warning(paste0("No features for '", names(geneSetList)[i], "' were found in 'rownames(inSCE)'. Excluding this gene set.", sep=""))
        }

        if(length(gs.diff) > 0 & length(gs.i) > 0) {
          warning(paste0("Some features were not found in 'rownames(inSCE)' for gene set '",
                         names(geneSetList)[i], "': ",
                         paste(gs.diff, collapse=","), sep=""))
        }

      } else {
        ## If IDs located in rowData
        ## Check for existence of location in rowData(inSCE)
        if(any(locations[i] == seq(ncol(rowData(inSCE))))) {
          temp.location <- as.numeric(locations[i])
        } else if (any(locations[i] == colnames(rowData(inSCE)))) {
          temp.location <- as.character(locations[i])
        } else {
          stop(paste0("'", locations[i], "' was not found in the column names of 'rowData(inSCE).'", sep=""))
        }

        gs.i <- intersect(geneSetList[[i]], rowData(inSCE)[,temp.location])
        gs.diff <- setdiff(geneSetList[[i]], rowData(inSCE)[,temp.location])

        if(length(gs.i) > 0) {
          gs.index <- which(rowData(inSCE)[,temp.location] %in% gs.i)
          geneSets[[names(geneSetList)[i]]] <- rownames(inSCE)[gs.index]
        } else {
          warning(paste0("No features for gene set '", names(geneSetList)[i],
                         "' were found in 'rowData(inSCE)' under column '",
                         temp.location, "'. Excluding this gene set.", sep=""))
        }

        if(length(gs.diff) > 0 & length(gs.i) > 0) {
          warning(paste0("Some features were not found in 'rowData(inSCE)' under column '",
                         temp.location, " 'for gene set '",
                         names(geneSetList)[i], "': ",
                         paste(gs.diff, collapse=","), sep=""))
        }
      }
    }
  }

  # ## Add mito genes into geneSets 
  # if (length(mitoG) != 0) {
  #   geneSets[['mito']] <- mitoG
  # }

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
    
  inSCE <- scater::addPerCellQC(x = inSCE,
                                exprs_values = useAssay,
                                subsets = geneSets,
                                percent_top = percent_top,
                                use_altexps = use_altexps,
                                flatten = flatten,
                                detection_limit = detectionLimit,
                                BPPARAM = BPPARAM)

  ## rename mito gene columns in colData(inSCE)
  mitoCols <- grep("subsets_mito", names(SummarizedExperiment::colData(inSCE)))
  newMitoCols <- gsub('subsets_mito', 'mito', names(SummarizedExperiment::colData(inSCE))[mitoCols])
  names(SummarizedExperiment::colData(inSCE))[mitoCols] <- newMitoCols

  argsList <- argsList[!names(argsList) %in% ("BPPARAM")]
  S4Vectors::metadata(inSCE)$scater_addPerCellQC <- argsList[-1]
  S4Vectors::metadata(inSCE)$scater_addPerCellQC$packageVersion <- utils::packageDescription("scran")$Version

  if(is.null(geneSets)){
    geneSets = as.character(geneSets)
  }
  S4Vectors::metadata(inSCE)$scater$addPerCellQC$geneSets <- geneSets

  return(inSCE)
}
