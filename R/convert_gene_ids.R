#' Convert gene IDs
#'
#' Convert the gene IDs in a SingleCellExperiment object using Bioconductor
#' org.*.eg.db data packages. Because annotation databases do not have a 1:1
#' relationship, this tool removes rows with no corresponding annotation in
#' your desired annotation, and remove any duplicate annotations after
#' conversion.
#'
#' @param inSCESet Input SCtkExperiment object. Required
#' @param in_symbol The input symbol type
#' @param out_symbol The output symbol type
#' @param database The org.*.eg.db database to use. The default is org.Hs.eg.db
#' 
#' @return A SCtkExperiment with converted gene IDs.
#' 
#' @export
convert_gene_ids <- function(inSCESet, in_symbol, out_symbol,
                             database="org.Hs.eg.db"){
  if(!(database %in% c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db",
    "org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db","org.Dr.eg.db",
    "org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db",
    "org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db",
    "org.Rn.eg.db","org.Sc.sgd.db","org.Ss.eg.db","org.Xl.eg.db"))){
    stop("The database you want to use, ", database, ", is not supported")
  }
  if(!(database %in% as.character(grep("^org\\.", installed.packages()[, "Package"], value=T)))){
    stop("The database you want to use, ", database, ", is not installed.")
  }
  if(!(database %in% (.packages()))){
    stop("You need to load the database to use it: library(", database, ")")
  }
  if(in_symbol == out_symbol){
    message("No conversion necessary.")
    return(inSCESet)
  }
  indb <- get(paste(database))
  oldids <- rownames(inSCESet)
  if(any(duplicated(oldids)) | any(is.na(oldids))){
    stop("problem with input IDs, duplicates or NAs found")
  }
  res <- AnnotationDbi::select(indb, keys=oldids, columns=c(out_symbol),
                               keytype = in_symbol)
  if(any(is.na(res[,out_symbol]))){
    message(sum(is.na(res[,out_symbol])), " ", out_symbol ,
            " are NA after conversion. Removing")
    resna <- res[is.na(res[,out_symbol]),]
    res <- res[!is.na(res[,out_symbol]),]
    if(any(is.na(res[,out_symbol]))){
      stop("problem.")
    }
  }
  if(any(duplicated(res[,out_symbol]))){
    message(sum(duplicated(res[,out_symbol])), " ", out_symbol,
            " are duplicated after conversion. Removing additional copies")
    resdup <- res[duplicated(res[,out_symbol]),]
    res <- res[!duplicated(res[,out_symbol]),]
    if(any(duplicated(res[,out_symbol]))){
      stop("problem.")
    }
  }
  if(any(duplicated(res[,in_symbol]))){
    message(sum(duplicated(res[,in_symbol])), " ", in_symbol,
            " are duplicated after conversion. Removing additional copies")
    resdup_orig <- res[duplicated(res[,in_symbol]),]
    res <- res[!duplicated(res[,in_symbol]),]
    if(any(duplicated(res[,in_symbol]))){
      stop("problem.")
    }
  }
  message(length(oldids), " ", in_symbol, " originally, " , nrow(res), " ",
          out_symbol, "s")
  newsce <- inSCESet[res[, in_symbol], ]
  rownames(newsce) <- res[, out_symbol]
  return(newsce)
}
