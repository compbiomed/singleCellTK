descriptionRunPerCellQC <- function() {
    return(list(
        introduction = "SingleCellTK utilizes the
    [scater](https://bioconductor.org/packages/release/bioc/html/scater.html)
    package to compute cell-level QC metrics.",
        runPerCellQC = "The wrapper function `runPerCellQC` can be used to separately
    compute QC metrics on its own.",
        parameter = "In this function, the `inSCE` parameter is the input
    SingleCellExperiment object, while the `useAssay` parameter is the assay
    object that in the SingleCellExperiment object the user wishes to use.",
        geneSet = "If the user wishes, a list of gene sets can be applied to the
    function to determine the expression of a set of specific genes.
    A gene list imported into the SingleCellExperiment
    object using `importGeneSets` functions can be set as `collectionName`.
    Additionally, a pre-made list of genes can be used to determine the level
    of gene expression per cell. A list containing gene symbols may be set as `geneSetList`,
    or the user may instead use the
    `geneSetCollection` parameter to supply a `GeneSetCollection` object from the
    [GSEABase](https://bioconductor.org/packages/release/bioc/html/GSEABase.html)
    package.",
        output = "The QC outputs are `sum`, `detected`, and `percent_top_X`.",
        sum = "`sum` contains the total number of counts for each cell.",
        detected = "`detected` contains the total number of features for each cell.",
        percentTop = "`percent_top_X` contains the percentage of the total counts that is made up
    by the expression of the top X genes for each cell.",
        subsets = "The `subsets_` columns contain information for the specific gene list that
    was used. For instance, if a gene list containing mitochondrial genes
    named `mito` was used, `subsets_mito_sum` would contains the total number of
    mitochondrial counts for each cell.",
        plotRunPerCellQCResults = "The wrapper function `plotRunPerCellQCResults` can be used
    to plot the general QC outputs."
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
        runEmptyDrops = "The wrapper function `runEmptyDrops` can be used to separately run the
            EmptyDrops algorithm on its own.",
        parameter = "In `runEmptyDrops`, the `lower` parameter is the lower bound of the
             total UMI count, in which all barcodes below the lower bound are
             assumed to be empty droplets. The `niters` parameter is the number
             of iterations the function will run for the calculation.
            `testAmbient` indicates whether results should be returned for
            barcodes that have a total UMI count below what is specified in
            `lower`.",
        plotEmptyDropsResults = "The wrapper function `plotEmptyDropsResults` can be used to plot the
              results from the EmptyDrops algorithm. This will visualize the empty droplets,
              by plotting the total UMI counts against the log probability for each barcode.",
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
        runBarcodeRankDrops = "The wrapper function `runBarcodeRankDrops` can be used to separately run the
              BarcodeRanks algorithm on its own.",
        parameter = "The `lower` parameter is again the lower bound of the total UMI
             count, in which all barcodes below the lower bound are assumed to
             be empty droplets.",
        plotBarcodeRankDropsResults = "The wrapper function `plotBarcodeRankDropsResults` can be used to plot the
              results from the BarcodeRanks algorithm.",
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
        runScrublet = "The wrapper function `runScrublet` can be used to separately run the
    Scrublet algorithm on its own.",
        plotScrubletResults = "The wrapper function `plotScrubletResults` can be used to plot the
    results from the Scrublet algorithm.",
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
        runDoubletFinder = "The wrapper function `runDoubletFinder` can be used to separately run the
            DoubletFinder algorithm on its own.",
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
            assignment of whether the cell is a doublet.",
        plotDoubletFinderResults = "The wrapper function `plotDoubletFinderResults` can be used to plot the
              QC outputs from the DoubletFinder algorithm."
    ))
}

descriptionDoubletCells <- function() {
    return(list(
        introduction = "[DoubletCells](https://rdrr.io/github/MarioniLab/scran/man/doubletCells.html) is a doublet detection algorithm in the `scran`
             package. DoubletCells aims to detect doublets by creating a
             simulated doublet from existing cells and projecting it to
             the same PCA space as the cells.",
        runDoubletCells = "The wrapper function `runBarcodeRankDrops` can be used to separately run the
              DoubletCells algorithm on its own.",
        parameter = "The `nNeighbors` parameter is the number of nearest neighbor
	used to calculate the density for doublet detection. `simDoublets` is used
	to determine the number of simulated doublets used for doublet detection.",
        output = "The output of DoubletCells is a `scran_doubletCells_score`.
             The doublet score of a droplet will be higher if the
             it is deemed likely to be a doublet.",
        plotDoubletCellsResults = "The wrapper function `plotDoubletCellsResults` can be used to plot the
              QC outputs from the DoubletCells algorithm."
    ))
}

descriptionCXDS <- function() {
    return(list(
        introduction = "CXDS, or co-expression based doublet scoring, is an algorithm in
             the [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package which employs a binomial model for the co-expression of
             pairs of genes to determine doublets.",
        runCxds = "The wrapper function `runCxds` can be used to separately run the
            CXDS algorithm on its own.",
        nTop = "In runCxds, the `ntop` parameter is the number of top variance
             genes to consider.",
        binThresh = "The `binThresh` parameter is the minimum counts
             a gene needs to have to be included in the analysis.",
        verb = "`verb` determines whether progress messages will be displayed
             or not.",
        retRes = "`retRes` will determine whether the gene pair results
             should be returned or not.",
        estNdbl = "The user may set the estimated number of doublets with `estNdbl`.",
        output = "The output of runCxds is the doublet score, `scds_cxds_score`.",
        plotCxdsResults = "The wrapper function `plotCxdsResults` can be used to plot the
              QC outputs from the CXDS algorithm."
    ))
}

descriptionBCDS <- function() {
    return(list(
        introduction = "BCDS, or binary classification based doublet scoring, is an
            algorithm in the [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package which uses a binary classification approach to determine doublets.",
        runBcds = "The wrapper function `runBcds` can be used to separately run the
            BCDS algorithm on its own.",
        nTop = "In runBcds, the `ntop` parameter is the number of top variance
             genes to consider.",
        srat = "The `srat` parameter is the ratio between
            original number of cells and simulated doublets.",
        nmax = "The `nmax` parameter is the maximum number of cycles that the algorithm
            should run through. If set to `tune`, this will be automatic.",
        varImp = "The `varImp` parameter determines if the variable importance
             should be returned or not.",
        output = "The output of runBcds is `scds_bcds_score`, which is the
             likelihood that a cell is a doublet.",
        plotBcdsResults = "The wrapper function `plotBCDSResults` can be used to plot the
              QC outputs from the BCDS algorithm."
    ))
}

descriptionScdsHybrid <- function() {
    return(list(
        introduction = "The CXDS-BCDS hybrid algorithm, uses both CXDS and BCDS
             algorithms from the
            [SCDS](https://www.bioconductor.org/packages/release/bioc/html/scds.html)
            package.",
        runCxdsBcdsHybrid = "The wrapper function `runCxdsBcdsHybrid` can be used to separately run the
            CXDS-BCDS hybrid algorithm on its own.",
        parameters = "All parameters from the `runBCDS` and `runBCDS` functions
             may be applied to this function in the `cxdsArgs` and `bcdsArgs`
             parameters, respectively.",
        output = "The output of runCxdsBcdsHybrid is the doublet score,
             `scds_hybrid_score`.",
        plotScdsHybridResults = "The wrapper function `plotScdsHybridResults` can be used to plot the
              QC outputs from the CXDS-BCDS hybrid algorithm."
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
        runDecontX = "The wrapper function `runDecontX` can be used to separately run the
            DecontX algorithm on its own.",
        output = "The outputs of `runDecontX` are `decontX_contamination` and
             `decontX_clusters`.",
        contamination = "`decontX_contamination` is a numeric vector which characterizes
             the level of contamination in each cell.",
        clustering = "Clustering is performed as part of the `runDecontX` algorithm.
             `decontX_clusters` is the resulting cluster assignment,
             which can also be labeled on the plot.",
        plotDecontXResults = "The wrapper function `plotDecontXResults` can be used to plot the
              QC outputs from the DecontX algorithm."
    ))
}

descriptionRunCellQC <- function() {
    return(list(
        introduction = "All of the droplet-based QC algorithms are able to be run under the wrapper
             function `runCellQC`. By default all possible QC algorithms will be run.",
        algorithms = "If users choose to only run a specific set of algorithms,
            they can specify which to run with the `algorithms` parameter."
    ))
}

descriptionRunDropletQC <- function() {
    return(list(
        introduction = "All droplet-based QC functions are able to be run under
            the wrapper function `runDropletQC`. By default all possible QC algorithms will be run.",
        algorithms = "If users choose to only run a specific set of algorithms,
            they can specify which to run with the `algorithms` parameter."
    ))
}

description_subsetSCECols <- function() {
    return(list(
        introduction = "SingleCellExperiment objects can be subset by its colData using
    `subsetSCECols`.",
        colData = 'The `colData` parameter takes in an expression in character vector form
    which will be used to identify a subset of columns using variables found in the
    colData of the SingleCellExperiment object. For example, if x is a numeric vector
    in colData, then "x < 5" will return all columns with x less than 5.',
        params = "The `index` parameter takes in a vector of indices which should be kept,
    while `bool` takes in a vector of TRUE or FALSE which should be the same length as
    the number of columns in the SingleCellExperiment object."
    ))
}

