import scrublet as scr
import scipy.io
import numpy as np
import os
import sys
import matplotlib
import subprocess
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging 
import inspect



plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

def run_Scrublet(counts_matrix,genes,output,expected_doublet_rate=0.06,min_counts=2,min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30,
	umap_n_neighbors=10, min_dist=0.3, angle=0.9, tsne_n_neighbors=5, n_iter=1000, umap_order_points=True, tsne_order_points=True, fa_order_points=True): 
	
	try:
		os.makedirs(output)
	except  FileExistsError:
		sys.stdout.write("Directory already exists. Overwriting previous scrublet output if any.\n")
		#sys.exit()
	
	## Create log file to track scrublet std out/err
	logging.basicConfig(filename=os.path.join(output,"run.log"),level=logging.INFO)
	
	kwarg_values = [expected_doublet_rate,min_counts,min_cells, min_gene_variability_pctl, n_prin_comps,umap_n_neighbors, min_dist, angle, tsne_n_neighbors, n_iter, umap_order_points, tsne_order_points, fa_order_points]
	kwargs = ["expected_doublet_rate","min_counts","min_cells", "min_gene_variability_pctl", "n_prin_comps","umap_n_neighbors", "min_dist", "angle", "tsne_n_neighbors", "n_iter", "umap_order_points", "tsne_order_points", "fa_order_points"]

	logging.info('Counts matrix provided: {}'.format(counts_matrix))
	logging.info('Genes provided: {}'.format(genes))
	logging.info('Scrublet Parameters provided: {}'.format([keys+":"+str(values) for keys , values in zip(kwargs,kwarg_values)]))

	counts_matrix = scipy.io.mmread(counts_matrix).T.tocsc()
	genes = np.array(scr.load_genes(genes, delimiter='\t', column=1))

	logging.info('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
	logging.info('Number of genes in gene list: {}'.format(len(genes)))
	logging.info('Writing all ouputs to folder:  {}'.format(output))
	
	logging.info('Creating scrublet histogram...')
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts,min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps)

	scrubhist_fig, scrubhist_ax  =  scrub.plot_histogram();
	scrubhist_fig.savefig(os.path.join(output,"scrublet_histogram.png"))

	logging.info('Running UMAP...')
	scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, n_neighbors=umap_n_neighbors, min_dist=min_dist))
	umapfig, umapax = scrub.plot_embedding('UMAP', order_points=umap_order_points);
	umapfig.savefig(os.path.join(output,"scrublet_umap.png"))

	logging.info('Running t-SNE...')
	scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=angle))
	tsnefig, tsneax = scrub.plot_embedding('tSNE', order_points=tsne_order_points);
	tsnefig.savefig(os.path.join(output,"scrublet_tsne.png"))

	logging.info('Running ForceAtlas2..')
	scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=tsne_n_neighbors, n_iter=n_iter))

	fafig, faax = scrub.plot_embedding('FA', order_points=fa_order_points);
	fafig.savefig(os.path.join(output,"scrublet_fa.png"))
	
''' Test ''' 
#mtx = "/restricted/projectnb/camplab/projects/2019-06-20_CifuentesD/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/matrix.mtx"
#genes = "/restricted/projectnb/camplab/projects/2019-06-20_CifuentesD/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/genes.tsv"
#outpath = "/restricted/projectnb/camplab/projects/2019_ChallengeProject/sbandyadka/test"
#run_Scrublet(mtx,genes,outpath,expected_doublet_rate=0.1,fa_order_points=False)
