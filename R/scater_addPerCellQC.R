#' @title Wrapper for calculating QC metrics with scater.
#' @description A wrapper function for \link[celda]{decontX}. Identify
#'  potential contamination from experimental factors such as ambient RNA.
#' @param inSCE Input \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param useAssay  A string specifying which assay in the SCE to use. Default 'counts'.
#' @param geneSetList List of gene sets to be quantified. The genes in the assays will be matched to the genes in the list based on \code{geneSetListLocation}.
#' @param geneSetListLocation Character or numeric vector. If set to 'rownames', then the genes in 'geneSetList' will be looked up in \code{rownames(inSCE)}.
#' If another character is supplied, then genes will be looked up in the column names of \code{rowData(inSCE)}. A character vector with the same length as \code{geneSetList} can be supplied if the IDs for different
#' gene sets are found in different places, including a mixture of 'rownames' and \code{rowData(inSCE)}. An integer or integer vector can be supplied to denote the column index in \code{rowData(inSCE)}. Default 'rownames'.
#' @param geneSetCollection Class of \code{GeneSetCollection} from package \code{GSEAbase}. The location of the gene IDs in \code{inSCE} should be in the \code{description} slot of each gene set and should follow the
#' same notation as \code{geneSetListLocation}. The function \link[GSEABase]{getGmt} can be used to read in gene sets from a GMT file. If reading a GMT file, the second column for each gene set should be the description denoting the location
#' of the gene IDs in \code{inSCE}. These gene sets will be included with those from \code{geneSetList} if both parameters are provided.
#' @param ... Additional arguments to pass to \link[scater]{addPerCellQC}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  cell QC metrics added to the \link[SummarizedExperiment]{colData} slot. If \code{geneSetList} or \code{geneSetCollection} are provided, then the rownames for each gene set will be saved in \code{metadata(inSCE)$scater$addPerCellQC$geneSets}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' mito.ix = grep("^MT-", rowData(sce)$feature_name)
#' geneSet <- list("Mito"=rownames(sce)[mito.ix])
#' sce <- runPerCellQC(sce, geneSetList = geneSet)
#' @export
runPerCellQC <- function(inSCE,
    useAssay = "counts",
    geneSetList = NULL,
    geneSetListLocation = "rownames",
    geneSetCollection = NULL,
    ...
) {

  message(paste0(date(), " ... Running 'perCellQCMetrics'"))
  #argsList <- as.list(formals(fun = sys.function(sys.parent()), envir = parent.frame()))
  argsList <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Add gene sets in 'geneSetCollection' to 'geneSetList', if available
  original.length <- length(geneSetList)
  geneSetCollectionLocation <- c()
  if(!is.null(geneSetCollection)) {

    ## Get the location where the gene set Ids are stored in SCE object
    geneSetCollectionLocation <- sapply(geneSetCollection, GSEABase::description)

    ## If blank/null/NA, then set to rownames by default
    ix <- geneSetCollectionLocation == "" || is.na(geneSetCollectionLocation) || is.null(geneSetCollectionLocation)
    geneSetCollectionLocation[ix] <- "rownames"

    ## Add gene sets to geneSetList
    geneSetList <- c(geneSetList, GSEABase::geneIds(geneSetCollection))
  }

  ## Set up locations for all lists
  if(length(geneSetListLocation) == 1) {
    locations <- c(rep(geneSetListLocation, original.length), geneSetCollectionLocation)
  } else if(length(geneSetListLocation) == original.length) {
    locations <- c(geneSetListLocation, geneSetCollectionLocation)
  } else {
    stop("'geneSetListLocation' needs to be of length 1 and will be used for all gene sets in 'geneSetList' or needs to be the length of 'geneSetList'.")
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
		          names(geneSetList)[i], "': ", paste(gs.diff, collapse=","), sep=""))
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
		  warning(paste0("No features for gene set '", names(geneSetList)[i], "' were found in 'rowData(inSCE)' under column '",
		          temp.location, "'. Excluding this gene set.", sep=""))
		}

		if(length(gs.diff) > 0 & length(gs.i) > 0) {
		  warning(paste0("Some features were not found in 'rowData(inSCE)' under column '",
		   temp.location, "'for gene set '", names(geneSetList)[i], "': ", paste(gs.diff, collapse=","), sep=""))
		}
	  }
	}
  }


  if(length(geneSets) == 0) {
    geneSets <- NULL
  }
  inSCE <- scater::addPerCellQC(x = inSCE, exprs_values = useAssay, subsets = geneSets, ...)

  argsList <- argsList[!names(argsList) %in% ("...")]
  dotList <- list(...)
  dotList <- dotList[!names(dotList) %in% c("BPPARAM")]
  argsList <- c(argsList, dotList)
  S4Vectors::metadata(inSCE)$scater$addPerCellQC <- argsList[-1]
  S4Vectors::metadata(inSCE)$scater$addPerCellQC$packageVersion <- utils::packageDescription("scran")$Version

  if(is.null(geneSets)){
      geneSets = as.character(geneSets)
  }
  S4Vectors::metadata(inSCE)$scater$addPerCellQC$geneSets <- geneSets

  return(inSCE)
}

