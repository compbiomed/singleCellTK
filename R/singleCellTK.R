#' Run the single cell analysis app
#'
#' Use this function to run the single cell analysis app.
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param includeVersion Include the version number in the SCTK header. The
#' default is TRUE.
#' @param theme The bootswatch theme to use for the singleCellTK UI. The default
#' is 'flatly'.
#'
#' @import GSVAdata Biobase DelayedArray
#'
#' @return The shiny app will open
#' @export
#' @examples
#' #Upload data through the app
#' if(interactive()){
#'   singleCellTK()
#' }
#'
#' # Load the app with a SingleCellExperiment object
#' if(interactive()){
#'   data("mouseBrainSubsetSCE")
#'   singleCellTK(mouseBrainSubsetSCE)
#' }
#'
singleCellTK <- function(inSCE=NULL, includeVersion=TRUE, theme='yeti') {
  appDir <- system.file("shiny", package = "singleCellTK")
  if (!is.null(inSCE) & is.null(rownames(inSCE))){
    stop("ERROR: No row names (gene names) found.")
  }
  shiny::shinyOptions(inputSCEset = inSCE)
  shiny::shinyOptions(includeVersion = includeVersion)
  shiny::shinyOptions(theme = theme)
  shiny::runApp(appDir, display.mode = "normal")
}
