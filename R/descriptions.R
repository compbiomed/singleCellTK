descriptionRunPerCellQC <- function() {
    return(list(
        introduction = "SingleCellTK utilizes the
    [scater](https://bioconductor.org/packages/release/bioc/html/scater.html)
    package to compute cell-level QC metrics in `runPerCellQC`.",
        parameter = "In this function, the `inSCE` parameter is the input
    SingleCellExperiment object, while the `useAssay` parameter is the assay
    object that in the SingleCellExperiment object the user wishes to use.",
        geneSet = "If the user wishes, a list of gene sets can be applied to the
    function to determine the expression of a set of specific genes in the
    `geneSetList` parameter. For instance, a pre-made list of mitochondrial
    genes can be used to determine the level of mitochondrial gene expression
    per cell. In lieu of the `geneSetList`, the user may instead use the
    `geneSetCollection` parameter to supply a `GeneSetCollection` object from the
    [GSEABase](https://bioconductor.org/packages/release/bioc/html/GSEABase.html)
    package.",
        output = "`runPerCellQC` produces the outputs `sum`, `detected`, and
    `percent_top_X`.",
        sum = "`sum` contains the total number of counts for each cell.",
        detected = "`detected` contains the total number of features for each cell.",
        percentTop = "`percent_top_X` contains the percentage of the total counts that is made up
    by the expression of the top X genes for each cell.",
        subsets = "The `subsets_` columns contain information for the specific gene list that
    was used. For instance, if a gene list containing mitochondrial genes
    named `mito` was used, `subsets_mito_sum` would contains the total number of
    mitochondrial counts for each cell."
    ))
}

descriptionEmptyDrops <- function() {
    return(list(
        introduction = "It is crucial to distinguish the data occurring from real cells
             and empty droplets containing ambient RNA. SCTK employs the
             [EmptyDrops](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html)
             algorithm from the
            [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
            package to test for empty droplets.",
        parameter = "In `runEmptyDrops`, the `lower` parameter is the lower bound of the
             total UMI count, in which all barcodes below the lower bound are
             assumed to be empty droplets. The `niters` parameter is the number
             of iterations the function will run for the calculation.
            `testAmbient` indicates whether results should be returned for
            barcodes that have a total UMI count below what is specified in
            `lower`.",
        plot = "To visualize the empty droplet, we will plot the total
             UMI counts against the log probability for each barcode.",
        plot2 = "Data points are colored by FDR values, where we see a small portion
            of the dataset contains barcodes that do not meet the threshold."
    ))
}

descriptionBarcodeRank <- function() {
    return(list(
        introduction = "[BarcodeRanks](https://rdrr.io/bioc/DropletUtils/man/barcodeRanks.html)
    from the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
            package computes barcode rank statistics and identifies
             the knee and inflection points on the total count curve. The knee
             and inflection points on the curve represent the difference between
             empty droplets and cell-containing droplets with much more RNA.",
        parameter = "The `lower` parameter is again the lower bound of the total UMI
             count, in which all barcodes below the lower bound are assumed to
             be empty droplets.",
        plot = "The total UMI count of each barcode is plotted against its rank, where
              we see a steep dropoff of UMI counts around the inflection point,
              where we see a separation between cell containing and empty droplets."
    ))
}

descriptionScrublet <- function() {
    return(list(
        introduction = "[Scrublet](https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb) aims to detect doublets by
             creating simulated doublets from combining transcriptomic profiles of existing cells in the dataset.",
        parameter = "The `sample` parameter indicates what sample each cell originated from.
             It can be set to `NULL` if all cells in the dataset came from the same sample.",
        additionalParam = "Scrublet also has a large set of parameters that the user can adjust;
             please refer to the Scrublet website for more details.",
        output = "The Scrublet outputs are `scrublet_score`, which is a numeric
            variable of the likelihood that a cell is a doublet, and the
            `scrublet_label`, which is the assignment of whether the cell
            is a doublet."
    ))
}

descriptionDoubletFinder <- function() {
    return(list(
        introduction = "[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) is a doublet detection algorithm which depends on
            the single cell analysis package
            [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html).",
        seuratRes = "`runDoubletFinder` relies on a parameter (in Seurat) called
                resolution to determine cells that may be doublets. Users
                will be able to manipulate the resolution parameter through
                `seuratRes`. If multiple numeric vectors are stored in
                `seuratRes`, there will be multiple label/scores.",
        seuratNfeatures = "The `seuratNfeatures` parameter determines the number of features
            that is used in the `FindVariableFeatures` function in Seurat.",
        seuratPCs = "`seuratPcs` parameter determines the number of dimensions used in
            the `FindNeighbors` function in Seurat.",
        seuratFormationRate = "The `formationRate`
            parameter is the estimated doublet detection rate in the dataset.
            aims to detect doublets by creating simulated doublets from
            combining transcriptomic profiles of existing cells in the dataset.",
        output = "The DoubletFinder outputs are `doubletFinder_doublet_score`,
            which is a numeric variable of the likelihood that a cell is a
            doublet, and the `doubletFinder_doublet_label`, which is the
            assignment of whether the cell is a doublet."
    ))
}

descriptionDoubletCells <- function() {
    return(list(
        introduction = "[DoubletCells](https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html) is a doublet detection algorithm in the `scran`
             package. DoubletCells aims to detect doublets by creating a
             simulated doublet from existing cells and projecting it to
             the same PCA space as the cells.",
        parameter = "The `nNeighbors` parameter is the number of nearest neighbor
	used to calculate the density for doublet detection. `simDoublets` is used
	to determine the number of simulated doublets used for doublet detection.",
        output = "The output of `runDoubletCells` is a `scran_doubletCells_Score`.
             The doublet score of a droplet will be higher if the
             it is deemed likely to be a doublet."
    ))
}

descriptionCXDS <- function() {
    return(list(
        introduction = "CXDS, or co-expression based doublet scoring, is an algorithm in
             the [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package which employs a binomial model for the co-expression of
             pairs of genes to determine doublets.",
        nTop = "In runCxds, the `ntop` parameter is the number of top variance
             genes to consider.",
        binThresh = "The `binThresh` parameter is the minimum counts
             a gene needs to have to be included in the analysis.",
        verb = "`verb` determines whether progress messages will be displayed
             or not.",
        retRes = "`retRes` will determine whether the gene pair results
             should be returned or not.",
        estNdbl = "The user may set the estimated number of doublets with `estNdbl`.",
        output = "The output of runCxds is the doublet score, `scds_cxds_score`."
    ))
}

descriptionBCDS <- function() {
    return(list(
        introduction = "BCDS, or binary classification based doublet scoring, is an
            algorithm in the [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package which uses a binary classification approach to determine doublets.",
        nTop = "In runBcds, the `ntop` parameter is the number of top variance
             genes to consider.",
        srat = "The `srat` parameter is the ratio between
            original number of cells and simulated doublets.",
        nmax = "The `nmax` parameter is the maximum number of cycles that the algorithm
            should run through. If set to `tune`, this will be automatic.",
        varImp = "The `varImp` parameter determines if the variable importance
             should be returned or not.",
        output = "The output of runBcds is `scds_bcds_score`, which is the
             likelihood that a cell is a doublet."
    ))
}

descriptionScdsHybrid <- function() {
    return(list(
        introduction = "runCxdsBcdsHybrid, uses both CXDS and BCDS
             algorithms from the
            [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package.",
        parameters = "All parameters from the `runBCDS` and `runBCDS` functions
             may be applied to this function in the `cxdsArgs` and `bcdsArgs`
             parameters, respectively.",
        output = "The output of runCxdsBcdsHybrid is the doublet score,
             `scds_hybrid_score`."
    ))
}

descriptionDecontX <- function() {
    return(list(
        introduction = "In droplet-based single cell technologies,
            ambient RNA that may have been released from apoptotic or
            damaged cells may get incorporated into another droplet, and can
            lead to contamination. [decontX](https://rdrr.io/bioc/celda/man/decontX.html),
            available from the [celda](https://bioconductor.org/packages/release/bioc/html/celda.html),
            is a Bayesian method for the identification of the contamination level at a cellular level.",
        output = "The outputs of `runDecontX` are `decontX_contamination` and
             `decontX_clusters`.",
        contamination = "`decontX_contamination` is a numeric vector which characterizes
             the level of contamination in each cell.",
        clustering = "Clustering is performed as part of the `runDecontX` algorithm.
             `decontX_clusters` is the resulting cluster assignment,
             which can also be labeled on the plot."
    ))
}

descriptionRunCellQC <- function() {
    return(list(
        introduction = "All of the above functions are able to be run under the wrapper
             function `runCellQC`. By default all of the functions will be run.",
        algorithms = "If users choose to only run a specific set of algorithms,
            they can specify which to run with the `algorithms` parameter."
    ))
}

descriptionRunDropletQC <- function() {
    return(list(
        introduction = "All of the above droplet-based functions are able to be run under
            the wrapper function `runDropletQC`. By default all of the functions will be run.",
        algorithms = "If users choose to only run a specific set of algorithms,
            they can specify which to run with the `algorithms` parameter."
    ))
}