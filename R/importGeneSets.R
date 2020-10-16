#' @title Imports gene sets from a GMT file
#' @description Converts a list of gene sets stored in a GMT file into a
#' \linkS4class{GeneSetCollection} and stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param file Character. Path to GMT file. See \link[GSEABase]{getGmt} for
#' more information on reading GMT files.
#' @param collectionName Character. Name of collection to add gene sets to.
#' If this collection already exists in \code{inSCE}, then these gene sets will
#' be added to that collection. Any gene sets within the collection with the
#' same name will be overwritten. Default \code{GeneSetCollection}.
#' @param by Character, character vector, or NULL. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{geneSetList} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. \code{by} can be a vector the same length as
#' the number of gene sets in the GMT file and the elements of the vector
#' can point to different locations within \code{inSCE}. Finally, \code{by}
#' can be \code{NULL}. In this case, the location of the gene identifiers
#' in \code{inSCE} should be saved in the description (2nd column)
#' of the GMT file. See \link{featureIndex} for more information.
#' Default \code{"rownames"}.
#' @param sep Character. Delimiter of the GMT file. Default \code{"\t"}.
#' @details The gene identifiers in gene sets in the GMT file will be
#' mapped to the rownames of \code{inSCE} using the \code{by} parameter and
#' stored in a \linkS4class{GeneSetCollection} object from package
#' \link{GSEABase}. This object is stored in
#' \code{metadata(inSCE)$sctk$genesets}, which can be accessed in downstream
#' analysis functions such as \link[singleCellTK]{runCellQC}.
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromList} for importing from lists,
#' \link{importGeneSetsFromCollection} for importing from
#' \linkS4class{GeneSetCollection} objects, and
#' \link{importGeneSetsFromMSigDB} for importing MSigDB gene sets.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  with gene set from \code{collectionName} output stored to the
#'  \link[S4Vectors]{metadata} slot.
#' @examples
#' data(scExample)
#'
#' # GMT file containing gene symbols for a subset of human mitochondrial genes
#' gmt <- system.file("extdata/mito_subset.gmt", package = "singleCellTK")
#'
#' # "feature_name" is the second column in the GMT file, so the ids will
#' # be mapped using this column in the 'rowData' of 'sce'. This
#' # could also be accomplished by setting by = "feature_name" in the
#' # function call.
#' sce <- importGeneSetsFromGMT(inSCE = sce, file = gmt, by = NULL)
#' @export
importGeneSetsFromGMT <- function(inSCE, file,
                                  collectionName = "GeneSetCollection",
                                  by = "rownames", sep = "\t") {

  # Read gmt file into GeneSetCollection
  gsc <- GSEABase::getGmt(con = file, sep = sep)

  inSCE <- importGeneSetsFromCollection(inSCE = inSCE,
                                        geneSetCollection = gsc,
                                        collectionName = collectionName,
                                        by = by)
  return(inSCE)
}


#' @title Imports gene sets from a list
#' @description Converts a list of gene sets into a
#' \linkS4class{GeneSetCollection} and stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param geneSetList Named List. A list containing one or more gene sets.
#' Each element of the list should be a character vector of gene identifiers.
#' The names of the list will be become the gene set names in the
#' \linkS4class{GeneSetCollection} object.
#' @param collectionName Character. Name of collection to add gene sets to.
#' If this collection already exists in \code{inSCE}, then these gene sets will
#' be added to that collection. Any gene sets within the collection with the
#' same name will be overwritten. Default \code{GeneSetCollection}.
#' @param by Character or character vector. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{geneSetList} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. Finally, \code{by} can be a vector the same length as
#' the number of gene sets in \code{geneSetList} and the elements of the vector
#' can point to different locations within \code{inSCE}. See
#' \link{featureIndex} for more information.
#' Default \code{"rownames"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object
#' with gene set from \code{collectionName} output stored to the
#' \link[S4Vectors]{metadata} slot.
#' @details The gene identifiers in gene sets in \code{geneSetList} will be
#' mapped to the rownames of \code{inSCE} using the \code{by} parameter and
#' stored in a \linkS4class{GeneSetCollection} object from package
#' \link{GSEABase}. This object is stored in
#' \code{metadata(inSCE)$sctk$genesets}, which can be accessed in downstream
#' analysis functions such as \link[singleCellTK]{runCellQC}.
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromCollection} for importing from
#' \linkS4class{GeneSetCollection} objects,
#' \link{importGeneSetsFromGMT} for importing from GMT files, and
#' \link{importGeneSetsFromMSigDB} for importing MSigDB gene sets.
#' @examples
#' data(scExample)
#'
#' # Generate gene sets from 'rownames'
#' gs1 <- rownames(sce)[1:10]
#' gs2 <- rownames(sce)[11:20]
#' gs <- list("geneset1" = gs1, "geneset2" = gs2)
#' sce <- importGeneSetsFromList(inSCE = sce,
#'                               geneSetList = gs,
#'                               by = "rownames")
#'
#' # Generate a gene set for mitochondrial genes using
#' # Gene Symbols stored in 'rowData'
#' mito.ix <- grep("^MT-", rowData(sce)$feature_name)
#' mito <- list(mito = rowData(sce)$feature_name[mito.ix])
#' sce <- importGeneSetsFromList(inSCE = sce,
#'                              geneSetList = mito,
#'                              by = "feature_name")
#' @export
#' @importFrom SummarizedExperiment rowData
importGeneSetsFromList <- function(inSCE, geneSetList,
                                   collectionName = "GeneSetCollection",
                                   by = "rownames") {

  # Check to ensure 'by' is correct
  if(!all(by %in% c("rownames", colnames(rowData(inSCE))))) {
    stop("The entries in 'by' must be 'rownames' or one of the column names ",
         "in the 'rowData(inSCE)': ", paste(colnames(rowData(inSCE)),
         collapse = ", "))
  }
  if(length(by) != 1 & length(by) != length(geneSetList)) {
    stop("'by' needs to be a character of length 1 describing the location ",
         "of all the gene sets in inSCE or a character vector the same ",
         "length as 'geneSetList' where each entry describes the location of ",
         "that particular gene set in 'geneSetList'.")
  }
  if(is.null(names(geneSetList))) {
    stop("'geneSetList' needs to be a named list.")
  }
  if(any(is.na(names(geneSetList)))) {
    stop("No names in 'geneSetList' can be 'NA'.")
  }

  # Convert to GeneSetCollection
  gs <- list()
  for(i in seq_along(geneSetList)) {
    gs[[i]] <- GSEABase::GeneSet(setName = names(geneSetList)[i],
                                 geneIds = geneSetList[[i]])
  }
  gsc <- GSEABase::GeneSetCollection(gs)

  inSCE <- importGeneSetsFromCollection(inSCE = inSCE,
                                        geneSetCollection = gsc,
                                        collectionName = collectionName,
                                        by = by)
  return(inSCE)
}

#' @title Imports gene sets from a GeneSetCollection object
#' @description Converts a list of gene sets stored in a
#' \linkS4class{GeneSetCollection} object and stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param geneSetCollection A \linkS4class{GeneSetCollection} object. See
#' \link[GSEABase]{GeneSetCollection} for more details.
#' @param collectionName Character. Name of collection to add gene sets to.
#' If this collection already exists in \code{inSCE}, then these gene sets will
#' be added to that collection. Any gene sets within the collection with the
#' same name will be overwritten. Default \code{GeneSetCollection}.
#' @param by Character, character vector, or NULL. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{geneSetCollection} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. \code{by} can be a vector the same length as
#' the number of gene sets in the \code{GeneSetCollection} and the elements of the vector
#' can point to different locations within \code{inSCE}. Finally, \code{by}
#' can be \code{NULL}. In this case, the location of the gene identifiers
#' in \code{inSCE} should be saved in the description slot for each gene set
#' in the \code{GeneSetCollection}.
#' See \link{featureIndex} for more information.
#' Default \code{"rownames"}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object
#' with gene set from \code{collectionName} output stored to the
#' \link[S4Vectors]{metadata} slot.
#' @details The gene identifiers in gene sets in the
#' \code{GeneSetCollection} will be mapped to the rownames of
#' \code{inSCE} using the \code{by} parameter and
#' stored in a \linkS4class{GeneSetCollection} object from package
#' \link{GSEABase}. This object is stored in
#' \code{metadata(inSCE)$sctk$genesets}, which can be accessed in downstream
#' analysis functions such as \link[singleCellTK]{runCellQC}.
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromList} for importing from lists,
#' \link{importGeneSetsFromGMT} for importing from GMT files, and
#' \link{importGeneSetsFromMSigDB} for importing MSigDB gene sets.
#' @examples
#' data(scExample)
#' library(GSEABase)
#' gs1 <- GeneSet(setName = "geneset1", geneIds = rownames(sce)[1:10])
#' gs2 <- GeneSet(setName = "geneset2", geneIds = rownames(sce)[11:20])
#' gsc <- GeneSetCollection(list(gs1, gs2))
#' sce <- importGeneSetsFromCollection(inSCE = sce,
#'                                     geneSetCollection = gsc,
#'                                     by = "rownames")
#' @export
importGeneSetsFromCollection <- function(inSCE, geneSetCollection,
                                         collectionName = "GeneSetCollection",
                                         by = "rownames") {

  # If 'by' is NULL, then the location will be derived from the description
  if(is.null(by)) {
    location <- unlist(lapply(1:length(geneSetCollection),
                   function(i) GSEABase::description(geneSetCollection[[i]])))
  } else {
    if(length(by) != 1 & length(by) != length(geneSetCollection)) {
      stop("If not NULL, then 'by' needs to be a character of ",
           "length 1 describing the location ",
           "of all the gene sets in inSCE or a character vector the same ",
           "length as 'geneSetList' where each entry describes the location of ",
           "that particular gene set in 'geneSetList'.")
    }

    if(length(by) == 1) {
      location <- rep(by, length(geneSetCollection))
    } else {
      location <- by
    }
  }

  # Check to ensure 'by' matches annotation in SCE
  if(!all(location %in% c("rownames", colnames(rowData(inSCE))))) {
    stop("The entries in 'by' must be 'rownames' or one of the column names ",
         "in the 'rowData(inSCE)': ", paste(colnames(rowData(inSCE)),
                                            collapse = ", "))
  }

  gs <- list()
  for(i in seq_along(geneSetCollection)) {

    temp.gs <- featureIndex(
      features = GSEABase::geneIds(geneSetCollection[[i]]),
      inSCE = inSCE,
      by = location[i],
      exactMatch = TRUE,
      removeNA = TRUE,
      errorOnNoMatch = FALSE,
      warningOnPartialMatch = FALSE)
    if(length(temp.gs) > 0) {
      gs.new <- geneSetCollection[[i]]
      GSEABase::geneIds(gs.new) <- rownames(inSCE)[temp.gs]
      gs <- c(gs, gs.new)
    }
  }

  if(length(gs) == 0) {
    stop("No gene sets were succesfully imported.")
  }

  # Add GeneSetCollection back to metadata
  new.gsc <- GSEABase::GeneSetCollection(gs)
  old.gsc <- NULL
  if(!is.null(S4Vectors::metadata(inSCE)$sctk$genesets)) {
    if(collectionName %in% names(S4Vectors::metadata(inSCE)$sctk$genesets)) {
      old.gsc <- S4Vectors::metadata(inSCE)$sctk$genesets[[collectionName]]
    }
  }

  if(!is.null(old.gsc)) {
    # Remove any old gene sets with same name as new gene sets
    old.gsc <- old.gsc[!(names(old.gsc) %in% names(new.gsc))]

    # Combine old and new GeneSetCollections and save back to metadata
    S4Vectors::metadata(inSCE)$sctk$genesets[[collectionName]] <-
      GSEABase::GeneSetCollection(c(old.gsc, new.gsc))
  } else {
    S4Vectors::metadata(inSCE)$sctk$genesets[[collectionName]] <- new.gsc
  }

  return(inSCE)
}



#' @title Imports gene sets from MSigDB
#' @description Gets a list of MSigDB gene sets stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param categoryIDs Character vector containing the MSigDB gene set ids.
#' The column \code{ID} in the table returned by \code{getMSigDBTable()} shows
#' the list of possible gene set IDs that can be obtained.
#' @param species Character. Species available can be found using the function
#' \code{\link[msigdbr]{msigdbr_show_species}}. Default \code{"Homo sapiens"}.
#' @param mapping Character. One of "gene_symbol", "human_gene_symbol", or
#' "entrez_gene". Gene identifiers to be used for MSigDB gene sets. IDs
#' denoted by the \code{by} parameter must be either in gene symbol or
#' Entrez gene id format to match IDs from MSigDB.
#' @param by Character. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' the MSigDB gene sets should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. See \link{featureIndex} for more information.
#' Default \code{"rownames"}.
#' @param verbose Boolean. Whether to display progress. Default \code{TRUE}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object
#' with gene set from \code{collectionName} output stored to the
#' \link[S4Vectors]{metadata} slot.
#' @details The gene identifiers in gene sets from MSigDB will be retrieved
#' using the \code{\link{msigdbr}} package. They will be mapped to the IDs in
#' \code{inSCE} using the \code{by} parameter and
#' stored in a \linkS4class{GeneSetCollection} object from package
#' \link{GSEABase}. This object is stored in
#' \code{metadata(inSCE)$sctk$genesets}, which can be accessed in downstream
#' analysis functions such as \link[singleCellTK]{runCellQC}.
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromList} for importing from lists,
#' \link{importGeneSetsFromGMT} for importing from GMT files, and
# \link{importGeneSetsFromCollection} for importing from
#' \linkS4class{GeneSetCollection} objects.
#' @examples
#' data(scExample)
#' sce <- importGeneSetsFromMSigDB(inSCE = sce,
#'                                 categoryIDs = "H",
#'                                 species = "Homo sapiens",
#'                                 mapping = "gene_symbol",
#'                                 by = "feature_name")
#' @export
#' @importFrom SummarizedExperiment rowData
importGeneSetsFromMSigDB <- function(inSCE, categoryIDs,
                                     species = "Homo sapiens",
                                     mapping = c("gene_symbol",
                                            "human_gene_symbol",
                                            "entrez_gene"),
                                     by = "rownames",
                                     verbose = TRUE) {

  mapping <- match.arg(mapping)

  # Check species
  all.species <- msigdbr::msigdbr_species()$species_name
  if(!(species %in% all.species)) {
    stop("'species' needs to be one of the following: ",
         paste(all.species, collapse = ", "))
  }
  # Check to ensure 'by' matches annotation in SCE
  if(is.null(by) || length(by) > 1) {
    stop("'by' must a character of length 1.")
  }
  if(!(by %in% c("rownames", colnames(rowData(inSCE))))) {
    stop("'by' must be 'rownames' or one of the column names ",
         "in the 'rowData(inSCE)': ", paste(colnames(rowData(inSCE)),
                                            collapse = ", "))
  }

  msigdb_table <- getMSigDBTable()
  if(!all(categoryIDs %in% msigdb_table$ID)) {
    diff <- setdiff(categoryIDs, msigdb_table$ID)
    stop("Category IDs not found in MSigDB: ", paste(diff, collapse = ", "),
         ". Valid IDs can be found in the table returned by the function ",
         "getMSigDBTable().")
  }

  rownames(msigdb_table) <- msigdb_table$ID
  gs.list <- list()
  for(i in seq_along(categoryIDs)) {

    category <- msigdb_table[categoryIDs[i],"Category"]
    subcat <- msigdb_table[categoryIDs[i],"Subcategory"]
    if(subcat == "N/A") {
      subcat <- NA
    }

    # Retrieve lists from msigdbr
    gs <- as.data.frame(msigdbr::msigdbr(species = species,
                    category = category,
                    subcategory = subcat))

    # Parse data.frame into list of GeneSets
    gs.names <- unique(gs$gs_name)
    num.gs.names <- length(gs.names)
    if(isTRUE(verbose)){
      message(paste0(date(), " .. Importing '", categoryIDs[i],
                     "' gene sets (n = ", num.gs.names, ")"))
    }

    for(j in seq_along(gs.names)) {
      gs.sub <- gs[gs$gs_name == gs.names[j],mapping]
      gs.list[[j]] <- GSEABase::GeneSet(setName = gs.names[j],
                        geneIds = gs.sub,
                        collectionType = GSEABase::BroadCollection(
                          category = tolower(category),
                          subCategory = subcat))
      if(j %% 1000 == 0) {
        if(isTRUE(verbose)) {
          message(paste0(date(), " .... Completed ", j, " out of ",
                         num.gs.names, " gene sets"))
        }
      }
    }
    if(isTRUE(verbose)) {
      message(paste0(date(), " .... Completed ", num.gs.names, " gene sets ",
                   "for ", i))
    }

    # Add to SCE
    if(isTRUE(verbose)) {
      message(paste0(date(), " .. Matching gene sets to '", by, "'"))
    }
    gsc <- GSEABase::GeneSetCollection(gs.list)
    inSCE <- importGeneSetsFromCollection(inSCE = inSCE,
                                          geneSetCollection = gsc,
                                          collectionName = categoryIDs[i],
                                          by = by)
  }

  return(inSCE)
}


#' @title Lists imported GeneSetCollections
#' @description Returns a vector of GeneSetCollections that have been
#' imported and stored in \code{metadata(inSCE)$sctk$genesets}.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @return Character vector.
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromList} for importing from lists,
#' \link{importGeneSetsFromGMT} for importing from GMT files,
# \link{importGeneSetsFromCollection} for importing from
#' \linkS4class{GeneSetCollection} objects, and \link{importGeneSetsFromMSigDB}
#' for importing MSigDB gene sets.
#' @examples
#' data(scExample)
#' library(GSEABase)
#' gs1 <- GeneSet(setName = "geneset1", geneIds = rownames(sce)[1:10])
#' gs2 <- GeneSet(setName = "geneset2", geneIds = rownames(sce)[11:20])
#' gsc1 <- GeneSetCollection(gs1)
#' gsc2 <- GeneSetCollection(gs2)
#' sce <- importGeneSetsFromCollection(inSCE = sce,
#'                                     geneSetCollection = gsc1,
#'                                     by = "rownames",
#'                                     collectionName = "Collection1")
#' sce <- importGeneSetsFromCollection(inSCE = sce,
#'                                     geneSetCollection = gsc2,
#'                                     by = "rownames",
#'                                     collectionName = "Collection2")
#' collections <- sctkListGeneSetCollections(sce)
#' @export
sctkListGeneSetCollections <- function(inSCE) {
  if(!is.null(S4Vectors::metadata(inSCE)$sctk$genesets)) {
    res <- names(S4Vectors::metadata(inSCE)$sctk$genesets)
  } else {
    res <- "No GeneSetCollections have been imported."
  }
  return(res)
}

#' @title Shows MSigDB categories
#' @description Returns a data.frame that shows MSigDB categories and
#' subcategories as well as descriptions for each. The entries in the ID
#' column in this table can be used as input for \link{importGeneSetsFromMSigDB}.
#' @return data.frame, containing MSigDB categories
#' @author Joshua D. Campbell
#' @seealso \link{importGeneSetsFromMSigDB} for importing MSigDB gene sets.
#' @export
getMSigDBTable <- function() {
  # Issues for getting this to pass R CMD Check:
  # https://support.bioconductor.org/p/24756/#24768
  # https://groups.google.com/forum/#!topic/cambridge-r-user-group/c7vf8o3QwDo
  msigdb_table <- NULL
  rm(msigdb_table)
  utils::data(msigdb_table)
  return(msigdb_table)
}


.retrieveGeneSetCollection <- function(inSCE, collectionName) {
  if(!is.null(S4Vectors::metadata(inSCE)$sctk$genesets)) {
    if(collectionName %in% names(S4Vectors::metadata(inSCE)$sctk$genesets)) {
      gsc <- S4Vectors::metadata(inSCE)$sctk$genesets[[collectionName]]
    } else {
      stop("No GeneSetCollection called '", collectionName, "' was found. ",
           "The GeneSetColletions that are available include: ",
           paste(names(S4Vectors::metadata(inSCE)$sctk$genesets),
                 collapse = ", "))
    }
  } else {
    stop("No GeneSetCollections have been imported.")
  }
  return(gsc)
}
