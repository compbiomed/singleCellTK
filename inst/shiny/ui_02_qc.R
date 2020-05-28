shinyPanelQC <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("Data QC"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),
    wellPanel(
      sidebarLayout(
        sidebarPanel(
          fluidRow(
            column(12, h3("Choose which algorithms to run:"))
          ),
          actionLink("selectallQC","Select All"),
          checkboxGroupInput("qcAlgos", "",
                             choiceNames =
                               list("doubletCells", "cxds", "bcds",
                                    "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder"),
                             choiceValues =
                               list("doubletCells", "cxds", "bcds",
                                    "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder")
          ),
          actionButton("runQC", "Run"),
          tags$div(id = "qcPageErrors"),
          # UNCOMMENT BELOW to display summary table after QC (must uncomment in server as well)
          # hidden(wellPanel(id = "qcData",
          #                  h4("Data summary:"),
          #                  tableOutput("qcSummary")))
        ),
        mainPanel(
        )
      )
    )
  )
)