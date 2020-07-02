#' Convert Gene IDs
#'
#' Convert the gene IDs in a SingleCellExperiment object using Bioconductor
#' org.*.eg.db data packages. Because annotation databases do not have a 1:1
#' relationship, this tool removes rows with no corresponding annotation in
#' your desired annotation, and remove any duplicate annotations after
#' conversion.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param inSymbol The input symbol type
#' @param outSymbol The output symbol type
#' @param database The org.*.eg.db database to use. The default is org.Hs.eg.db
#'
#' @return A \linkS4class{SingleCellExperiment} object with converted gene IDs.
#' @export
#' @examples
#' if(requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
#'   #convert mouse gene symbols to ensembl IDs
#'   library("org.Mm.eg.db")
#'   sample(rownames(mouseBrainSubsetSCE), 50)
#'   mouseBrainSubsetSymbol <- convertGeneIDs(inSCE = mouseBrainSubsetSCE,
#'                                            inSymbol = "SYMBOL",
#'                                            outSymbol = "ENSEMBL",
#'                                            database = "org.Mm.eg.db")
#'   sample(rownames(mouseBrainSubsetSymbol), 50)
#' }
#'
convertGeneIDs <- function(inSCE, inSymbol, outSymbol, database="org.Hs.eg.db"){
  if (!(database %in% c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db",
    "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
    "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db",
    "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db",
    "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db"))){
    stop("The database you want to use, ", database, ", is not supported")
  }
  if (!(database %in% as.character(grep("^org\\.",
                                      utils::installed.packages()[, "Package"],
                                       value = TRUE)))){
    stop("The database you want to use, ", database, ", is not installed.")
  }
  if (!(database %in% (.packages()))){
    stop("You need to load the database to use it: library(", database, ")")
  }
  if (inSymbol == outSymbol){
    message("No conversion necessary.")
    return(inSCE)
  }
  indb <- get(paste(database))
  oldids <- rownames(inSCE)
  if (any(duplicated(oldids)) | any(is.na(oldids))){
    stop("problem with input IDs, duplicates or NAs found")
  }
  res <- AnnotationDbi::select(indb, keys = oldids, columns = c(outSymbol),
                               keytype = inSymbol)
  if (any(is.na(res[, outSymbol]))){
    message(sum(is.na(res[, outSymbol])), " ", outSymbol,
            " are NA after conversion. Removing")
    res <- res[!is.na(res[, outSymbol]), ]
    if (any(is.na(res[, outSymbol]))){
      stop("problem.")
    }
  }
  if (any(duplicated(res[, outSymbol]))){
    message(sum(duplicated(res[, outSymbol])), " ", outSymbol,
            " are duplicated after conversion. Removing additional copies")
    res <- res[!duplicated(res[, outSymbol]), ]
    if (any(duplicated(res[, outSymbol]))){
      stop("problem.")
    }
  }
  if (any(duplicated(res[, inSymbol]))){
    message(sum(duplicated(res[, inSymbol])), " ", inSymbol,
            " are duplicated after conversion. Removing additional copies")
    res <- res[!duplicated(res[, inSymbol]), ]
    if (any(duplicated(res[, inSymbol]))){
      stop("problem.")
    }
  }
  message(length(oldids), " ", inSymbol, " originally, ", nrow(res), " ",
          outSymbol, "s")
  newsce <- inSCE[res[, inSymbol], ]
  rownames(newsce) <- res[, outSymbol]
  return(newsce)
}
