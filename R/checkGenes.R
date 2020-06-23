#' The checkGene function checks if submitted genes are recognized by a
#' database in EnrichR, then returns a list of either recognized or
#' unrecognized genes. 
#'
#' @param genelist A character vector of gene names to be uploaded to EnrichR.
#' @param database A character vector of reference database or databases to
#' match the genelist to.
#' @param baseurl Base url of EnrichR database gene libraries. Appended with
#' specific gene library name. Default:
#' "https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName="
#' @param output "recognized" to output the genes recognized in the gene set.
#' "unrecognized" to output the genes not recognized by the reference gene set.
#' Default: "recognized"
#' 
#' @return checkGene(): A list of recognized or unrecognized genes.
#' @export
#' @examples
#' databases=c("Reactome_2016", "KEGG_2019_Human", "WikiPathways_2019_Human", "Panther_2016")
#' genes=c("A1BG", "A1CF", "FAKEGENE")
#' checkGene(genes, databases, "recognized")
#' output= c("A1BG", "A1CF")
#' 
#' checkGene(genes, databases, "unrecognized")
#' output= c("FAKEGENE")

checkGene <- function(genelist, database, baseurl=
		      "https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName="
	      output="recognized") {
  ## Ensure databases are characters
  database=as.character(database)
 
  
  ## initialize list (to save time)
  database_loc <- list()
  database_txt <- list()
  no_col <- list()
  df <- list()
  
  ## import database data downloaded from Enrichr and organize into a df
  for(db in 1:length(database)) {
    database_loc[db] <- paste0(baseurl,database[db])
  
    ## To properly read the txt, read.table requires number of columns or
    ## variables (genes, pathways here are observations or rows)
    no_col[[db]] <- max(count.fields(database_loc[[db]], sep = "\t"), na.rm = TRUE)
    database_txt[[db]] <-
	    read.table(database_loc[[db]],sep="\t",fill=TRUE,header =
		       F,col.names=paste0("V",seq_len(no_col[[db]])),quote="")
  
    
    ## make a usable dataframe and list **might be a more efficient way**
    df[[db]] <- as.data.frame(database_txt[[db]][,-2])
    df[[db]] <- column_to_rownames(df[[db]], var = "V1")
    df[[db]] <- unite(df[[db]], col="Genes", sep=";")
  
    
    ## genelist
    df[[db]] <- sapply(df[[db]], strsplit, ";")
  
    ## get rid of gene duplicates
    df[[db]] <- unique(unlist(df[[db]]))
    
  }
  
  ## Character vector of all unique genes in databases
  path_lib <- unique(unlist(df))
  
  ## check if genes are in database list
  
  rec_genes <- genelist[genelist %in% path_lib]
  unrec_genes <- genelist[!genelist %in% path_lib]
  
  ##return either recognized genes or unrecognized genes
  if(output=="recognized"){
    return(rec_genes)
  }
  else if(output=="unrecognized") {
    return(unrec_genes)
  }
  else {
    print("Please define output as 'recognized' or 'unrecognized'. ")
  }
}

