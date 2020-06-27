#' @title Imports gene sets from a GMT file
#' @description Converts a list of gene sets stored in a GMT file into a
#' \linkS4class{GeneSetCollection} and stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param file Character. Path to GMT file. See \link[GSEABase]{getGmt} for
#' more information on reading GMT files. 
#' @param by Character, character vector, or NULL. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{GeneSetList} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. \code{'by'} can be a vector the same length as
#' the number of gene sets in the GMT file and the elements of the vector
#' can point to different locations within \code{inSCE}. Finally, \code{'by'}
#' can be \code {NULL}. In this case, the location of the gene identifiers 
#' in \code{inSCE} should be saved in the description (2nd column)
#' of the GMT file. See \link[celda]{retrieveFeatureIndex} for more information.
#' Default \code{"rownames"}.
#' @param sep Charcter. Delimiter of the GMT file. Default \code{"\t"}.
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
importGeneSetsFromGMT <- function(inSCE, file,
                                 by = "rownames", sep = "\t") {

  # Read gmt file into GeneSetCollection
  gsc <- GSEABase::getGmt(con = file, sep = sep)
  
  inSCE <- importGeneSetsFromCollection(inSCE = inSCE,
                                        geneSetCollection = gsc,
                                        by = by)
  return(inSCE)
}


#' @title Imports gene sets from a list
#' @description Converts a list of gene sets into a
#' \linkS4class{GeneSetCollection} and stores it in the metadata of the
#' \linkS4class{SingleCellExperiment} object. These gene sets can be used in
#' downstream quality control and analysis functions in \link{singleCellTK}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param GeneSetList Named List. A list containing one or more gene sets.
#' Each element of the list should be a character vector of gene identifiers. 
#' The names of the list will be become the gene set names in the 
#' \linkS4class{GeneSetCollection} object.
#' @param by Character or character vector. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{GeneSetList} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. Finally, \code{'by'} can be a vector the same length as
#' the number of gene sets in \code{GeneSetList} and the elements of the vector
#' can point to different locations within \code{inSCE}. See 
#' \link[celda]{retrieveFeatureIndex} for more information.
#' Default \code{"rownames"}.
#' @details The gene identifiers in gene sets in \code{GeneSetList} will be
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
#'                               GeneSetList = gs,
#'                               by = "rownames")
#' 
#' # Generate a gene set for mitochondrial genes using
#' # Gene Symbols stored in 'rowData'
#' mito.ix <- grep("^MT-", rowData(sce)$feature_name)
#' mito <- list(mito = rowData(sce)$feature_name[mito.ix])
#' sce <- importGeneSetsFromList(inSCE = sce,
#'                              GeneSetList = mito,
#'                              by = "feature_name")
importGeneSetsFromList <- function(inSCE, GeneSetList, by = "rownames") {

  # Check to ensure 'by' is correct
  if(!all(by %in% c("rownames", colnames(rowData(inSCE))))) {
    stop("The entries in 'by' must be 'rownames' or one of the column names ", 
         "in the 'rowData(inSCE)': ", paste(colnames(rowData(inSCE)), 
         collapse = ", "))
  }
  if(length(by) != 1 & length(by) != length(GeneSetList)) {
    stop("'by' needs to be a character of length 1 describing the location ",
         "of all the gene sets in inSCE or a character vector the same ",
         "length as 'GeneSetList' where each entry describes the location of ",
         "that particular gene set in 'GeneSetList'.")
  }
  if(is.null(names(GeneSetList))) {
    stop("'GeneSetList' needs to be a named list.")
  }
  if(any(is.na(names(GeneSetList)))) {
    stop("No names in 'GeneSetList' can be 'NA'.")
  }
  
  if(length(by) == 1) {
    location <- rep(by, length(GeneSetList))
  } else {
    location <- by
  }
  
  gs <- list()
  for(i in seq_along(GeneSetList)) {
    
    temp.gs <- suppressWarnings(celda::retrieveFeatureIndex(
      features = GeneSetList[[i]],
      x = inSCE,
      by = location[i],
      exactMatch = TRUE,
      removeNA = TRUE))
    if(length(temp.gs) > 0) {
      gs <- c(gs, GSEABase::GeneSet(setName = names(GeneSetList)[i],
                                    geneIds = rownames(inSCE)[temp.gs]))
    }
  }
  
  if(length(gs) == 0) {
    stop("No gene sets were succesfully imported.")
  }
  
  # Add GeneSetCollection back to metadata
  new.gsc <- GSEABase::GeneSetCollection(gs)
  old.gsc <- S4Vectors::metadata(inSCE)$sctk$genesets
  if(!is.null(old.gsc)) {
    # Remove any old gene sets with same name as new gene sets
    old.gsc <- old.gsc[!(names(old.gsc) %in% names(new.gsc))]
    
    # Combine old and new GeneSetCollections and save back to metadata
    S4Vectors::metadata(inSCE)$sctk$genesets <-
      GSEABase::GeneSetCollection(c(old.gsc, new.gsc))
  } else {
    S4Vectors::metadata(inSCE)$sctk$genesets <- new.gsc
  }
  
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
#' @param by Character, character vector, or NULL. Describes the
#' location within \code{inSCE} where the gene identifiers in
#' \code{GeneSetList} should be mapped. If set to \code{"rownames"} then the
#' features will be searched for among \code{rownames(inSCE)}. This can also be
#' set to one of the column names of \code{rowData(inSCE)} in which case the
#' gene identifies will be mapped to that column in the \code{rowData}
#' of \code{inSCE}. \code{'by'} can be a vector the same length as
#' the number of gene sets in the GMT file and the elements of the vector
#' can point to different locations within \code{inSCE}. Finally, \code{'by'}
#' can be \code {NULL}. In this case, the location of the gene identifiers 
#' in \code{inSCE} should be saved in the description slot for each gene set
#' in the \linkS4class{GeneSetCollection}.
#' See \link[celda]{retrieveFeatureIndex} for more information.
#' Default \code{"rownames"}.
#' @details The gene identifiers in gene sets in the
#' \linkS4class{GeneSetCollection will be mapped to the rownames of
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
importGeneSetsFromCollection <- function(inSCE, geneSetCollection,
                                         by = "rownames") {
  ids <- GSEABase::geneIds(geneSetCollection)
  
  # If 'by' is NULL, then the location will be derived from the description
  if(is.null(by)) {
    location <- unlist(lapply(1:length(geneSetCollection), 
                   function(i) GSEABase::description(geneSetCollection[[i]])))
  } else {
    location <- by
  }
  
  inSCE <- importGeneSetsFromList(inSCE = inSCE,
                                  GeneSetList = ids,
                                  by = location)
  return(inSCE)
}


# One of "C1" "C2" "C3" "C4"  "C5" "C6" "C7" "H"
importGeneSetsFromMSigDB <- function(inSCE, categoryIDs,
                                     mapping = c("gene_symbol",
                                            "human_gene_symbol",
                                            "entrez_gene"), 
                                     by = "rownames") {
  
  mapping <- match.arg(mapping)
  
  # Check species
  all.species <- msigdbr::msigdbr_show_species()
  if(!(species %in% all.species)) {
    stop("'species' needs to be one of the following: ",
         paste(all.species, collapse = ", "))
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
    for(j in gs.names) {
      gs.list[[j]] <- gs[gs$gs_name == j,mapping]
      #gs.sub <- gs[gs$gs_name == j,mapping]
      #gs.list[[j]] <- GSEABase::GeneSet(setName = j,
      #                  geneIds = gs.sub,
      #                  collectionType = GSEABase::BroadCollection(
      #                    category = tolower(category),
      #                    subCategory = subcat))
    }
  }  
  
  # Add to SCE
  gsc <- GSEABase::GeneSetCollection(gs.list)
  inSCE <- importGeneSetsFromCollection(inSCE = inSCE,
                                        geneSetCollection = gsc,
                                        by = by)
  return(inSCE)
}

#' @title Shows MSigDB categories
#' @description Returns a data.frame that shows MSigDB categories and 
#' subcategories as well as descriptions for each. The entries in the ID
#' column in this table can be used as input for \link{importGeneSetsFromMSigDB}.
#' @author Joshua D. Campbell
getMSigDBTable <- function() {
  data("msigdb_table")
  return(msigdb_table)
}