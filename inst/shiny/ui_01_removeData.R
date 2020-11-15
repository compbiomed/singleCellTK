shinyPanelRemove <- fluidPage(
  includeCSS('styles.CSS'),
  
  fluidRow(
    column(4,
           fluidRow(
             column(12,
                    panel(
                      heading = "Options",
                      selectInput(
                        inputId = "rmDataTypeSelect",
                        label = "Select type of data:",
                        choices = c("assays", "reducedDims"),
                        selected = "assays"
                      ),
                      selectInput("assayType", "Select assay:",
                                  list(
                                       "raw" = c("counts", "fpkm_counts"),
                                       "logCounts" = c("logCounts", "logSeuratCounts"),
                                       "normalized" = c("normCounts","seuratNormCounts"),
                                       "scaled" = c("seuratScaledCounts", "scaledCounts")
                                       )
                                  ),
                      withBusyIndicatorUI(
                        actionButton(
                          inputId = "delRedDim", 
                          label = "Delete")
                        )
                      )
                    )
           )
           ),
    column(8,
           fluidRow(
             column(12,
                    panel(
                      heading = "Data Available",
                      uiOutput(
                        outputId = "assaysList"
                      ),
                      uiOutput(
                        outputId = "reducedDimsList"
                        )
                      )
                    )
           )
           )
  )
)