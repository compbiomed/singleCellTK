#' @export
reportDifferentialFET <-
    function(sce,
             variable,
             phenotype,
             control,
             case,
             output_dir = ".",
             output_file = "DifferentialAbundanceFET_Report",
             pdf = FALSE,
             showSession = TRUE
             ){
    res <- differentialAbundanceFET(sce, variable, phenotype, control, case)
    rmarkdown::render(
        system.file("rmarkdown/DifferentialAbundanceFET_Report.Rmd"),
        params = list(
            sce = sce,
            res = res,
            pdf = isTRUE(pdf),
            showSession = isTRUE(showSession)
        ),
        output_file = output_file,
        output_dir = output_dir
    )
}


#' @export
reportClusterAbundance <- function(sce,
                                   variable,
                                   phenotype,
                                   output_dir = ".",
                                   output_file = "plotClusterAbundance_Report",
                                   pdf = FALSE,
                                   showSession = TRUE
                                   ){
    plot <- plotClusterAbundance(sce, variable, phenotype)
    rmarkdown::render(
        system.file("rmarkdown/PlotClusterAbundance_Report.Rmd"),
        params = list(
            sce = sce,
            plot = plot,
            pdf = isTRUE(pdf),
            showSession = isTRUE(showSession)
        ),
        output_file = output_file,
        output_dir = output_dir
    )
}

#' @export
reportMarkerDiffExp <- function(sce,
                                analysisName,
                                output_dir = ".",
                                output_file = "MarkerDiffExp_Report",
                                pdf = FALSE,
                                showSession = TRUE
                                ){
    de <- metadata(sce)$findMarker
    rmarkdown::render(
        system.file("rmarkdown/FindMarkerDiffExp_Report.Rmd"),
        params = list(
            sce = sce,
            de = de,
            output_file = output_file,
            pdf = isTRUE(pdf),
            showSession = isTRUE(showSession)
        ),
        output_file = output_file,
        output_dir = output_dir
    )
}