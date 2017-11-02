#' Run the single cell analysis app
#'
#' Use this function to run the single cell analysis app.
#'
#' @param inputData The input SCtkExperiment class object
#'
#' @return The shiny app will open
#'
#' @export singleCellTK
singleCellTK <- function(inputData=NULL) {
  appDir <- system.file("shiny", package = "singleCellTK")
  shiny::shinyOptions(inputSCEset = inputData)
  shiny::runApp(appDir, display.mode = "normal")
}
