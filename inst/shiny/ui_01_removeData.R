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
                      selectInput(
                        inputId = "delRedDimType", 
                        label = NULL, 
                        choices = currreddim),
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
                      heading = "Available:",
                      tableOutput(
                        outputId = "reducedDimsList"
                        )
                      )
                    )
           )
           )
  )
)