#' Run the single cell analysis app
#'
#' Use this function to run the single cell analysis app.
#'
#' @param inputData The input SCtkExperiment class object
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
#' #Load the app with a SCtkExperiment object
#' if(interactive()){
#'   data("mouseBrainSubsetSCE")
#'   singleCellTK(mouseBrainSubsetSCE)
#' }
#'
singleCellTK <- function(inputData=NULL) {
  appDir <- system.file("shiny", package = "singleCellTK")
  if (!is.null(inputData) & is.null(rownames(inputData))){
    stop("ERROR: No row names (gene names) found.")
  }
  shiny::shinyOptions(inputSCEset = inputData)
  shiny::runApp(appDir, display.mode = "normal")
}
