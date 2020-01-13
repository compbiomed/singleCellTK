#' @title Estimates doublet rates using Scrublet
#' @description Builds k-nearest neighbors using simulated doublets to estimate likelihood of doublets
#' @param counts_matrix full path to a counts matrix file - raw or filtered
#' @param genes full path to a gene list usually included in the raw sequecing folder
#' @param output full path to an output directory of preference. All plots, output files and log files will be written here
#' @param ... Any of expected_doublet_rate,min_counts,min_cells, min_gene_variability_pctl, n_prin_comps,umap_n_neighbors, min_dist, angle,
# tsne_n_neighbors, n_iter, umap_order_points, tsne_order_points, fa_order_points. Provided optionally as keyword arguments. Refer scrublet documentation for details
#' @return NA
#' @examples
#' runScrublet(counts_matrix=mtx,genes=genes,output=outpath,expected_doublet_rate=0.2)
#' @export
#' @import reticulate
### Requires python3/3.6.9 and  R/3.6.0

runScrublet <- function(counts_matrix,genes,output,...)
{	

	## Get user defined keyword arguments. If undefined, set defaults
	mc <- as.list(sys.call())
	if (!is.null(mc[['expected_doublet_rate']])) { edr <- mc$expected_doublet_rate } else { edr <- 0.06 }
	if (!is.null(mc[['min_counts']])) { mincounts <- mc$min_counts } else { mincounts <- 2 }
	if (!is.null(mc[['min_cells']])) { mincells <- mc$min_cells } else { mincells <- 3 }
	if (!is.null(mc[['min_gene_variability_pctl']])) { mingenevarpctl <- mc$min_gene_variability_pctl } else { mingenevarpctl <- 85 }
	if (!is.null(mc[['n_prin_comps']])) { nprincomps <- as.integer(mc$n_prin_comps) } else { nprincomps <- as.integer(30) }
	if (!is.null(mc[['umap_n_neighbors']])) { umapneighbors <- as.integer(mc$umap_n_neighbors) } else { umapneighbors <- as.integer(10) }
	if (!is.null(mc[['min_dist']])) { mindist <- mc$min_dist } else { mindist <- 0.3 }
	if (!is.null(mc[['angle']])) { ang <- mc$angle } else { ang <- 0.9 }
	if (!is.null(mc[['tsne_n_neighbors']])) { tsneneighbors <- as.integer(mc$tsne_n_neighbors) } else { tsneneighbors <- as.integer(5) }
	if (!is.null(mc[['n_iter']])) { niter <- mc$n_iter } else { niter <- 1000 }
	if (!is.null(mc[['umap_order_points']])) { umaporderpoints <- mc$umap_order_points } else { umaporderpoints <- 'True' }
	if (!is.null(mc[['tsne_order_points']])) { tsneorderpoints <- mc$tsne_order_points } else { tsneorderpoints <- 'True' }
	if (!is.null(mc[['fa_order_points']])) { faorderpoints <- mc$fa_order_points } else { faorderpoints <- 'True' }
	

	# Install reticulate library if not present.
	if("reticulate" %in% rownames(installed.packages()) == FALSE) {install.packages("reticulate")}
	library(reticulate)
	
	# check if required python packages in appropriate virtual env already exists. Create one if it doesn't 
	all_venvs <- virtualenv_list()
	if (grepl(pattern="r-scrubletvenv",x=all_venvs)==FALSE){
		virtualenv_create("r-scrubletvenv")
		virtualenv_install("r-scrubletvenv", "scipy")
		virtualenv_install("r-scrubletvenv", "scrublet")
		virtualenv_install("r-scrubletvenv", "matplotlib")
	}
	
	# Load virtual env and python script
	reticulate::use_virtualenv("r-scrubletvenv")
	reticulate::source_python("runscrublet.py")
	
	# call python scrublet function and write detected doublet scores to file
	
	doubletscore <- py_capture_output(run_Scrublet(counts_matrix,genes,output,expected_doublet_rate=edr,min_counts=mincounts,min_cells=mincells, min_gene_variability_pctl=mingenevarpctl, n_prin_comps=nprincomps,
	umap_n_neighbors=umapneighbors, min_dist=mindist, angle=ang, tsne_n_neighbors=tsneneighbors, n_iter=niter, umap_order_points=umaporderpoints, tsne_order_points=tsneorderpoints, fa_order_points=faorderpoints),
	 type = c("stderr","stdout"))
	
	doubletscore <- strsplit(doubletscore,'\n')[[1]]
	for (line in doubletscore){ 
		if (grepl("threshold",line)==TRUE) {doublet_score_threshold <- strsplit(line,"=")[[1]][2]}
		if (grepl("Detected",line)==TRUE) {detected_doublet_rate <- strsplit(line,"=")[[1]][2]}
		if (grepl("fraction",line)==TRUE) {detected_doublet_fraction <- strsplit(line,"=")[[1]][2]}
		if (grepl("\tExpected",line)==TRUE) {overall_expected_doublet_rate <- strsplit(line,"=")[[1]][2]}
		if (grepl("\tEstimated",line)==TRUE) {overall_estimated_doublet_rate <- strsplit(line,"=")[[1]][2]}
		if (grepl("Elapsed",line)==TRUE) {elapsed_time <- strsplit(line,":")[[1]][2]}
	}

	scrublet_estimates <- c(doublet_score_threshold,detected_doublet_rate,detected_doublet_fraction,overall_expected_doublet_rate,overall_estimated_doublet_rate,elapsed_time)
	
	doubletscore_outputfile <- file.path(output,"doublet_scores.txt")
	write(doubletscore, file=doubletscore_outputfile)

	return(scrublet_estimates)

}

#### Test #### 

#scrublet_source = "/restricted/projectnb/camplab/projects/2019_ChallengeProject/sbandyadka/"
#mtx = "/restricted/projectnb/camplab/projects/2019-06-20_CifuentesD/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/matrix.mtx"
#genes = "/restricted/projectnb/camplab/projects/2019-06-20_CifuentesD/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/genes.tsv"
#outpath = "/restricted/projectnb/camplab/projects/2019_ChallengeProject/sbandyadka/test"
#runScrublet(counts_matrix=mtx,genes=genes,output=outpath,expected_doublet_rate=0.06)
