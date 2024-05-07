renameClusterUI <- function(id) {
  ns <- NS(id)
  #print(table(old_label, new_label))
  #edit cluster labels in the dropdown
  #look at updated and previous and see the difference
  column(4, dropdown(
    fluidRow(
      column(width = 12,
             fluidRow(actionBttn(inputId = ns("closeDropDownFS1"), 
                                 label = NULL, style = "simple", 
                                 color = "danger",
                                 icon = icon("times"), size = "xs"), 
                      align = "right"),
             selectInput(
               inputId = ns("ClusterLabelChoices"),
               label = "Cluster Label",
               choices = NULL
             ),
             uiOutput(ns("textFieldsContainer")),
             textInput(
               inputId = ns("NewClusterLabel"),
               label = "Create a new cluster label (optional):",
             ),
             textOutput(ns("outputText")),
             withBusyIndicatorUI(
               actionBttn(
                 inputId = ns("clustPlot1"),
                 label = "Update",
                 style = "bordered",
                 color = "primary",
                 size = "sm"
               )
             ),
      )
    ),
    inputId = ns("dropDownFS1"),
    icon = icon("pencil"),
    status = "primary",
    circle = FALSE,
    inline = TRUE
  ))
}

renameClusterServer <- function(input, output, session, inSCE = NULL, isCluster = FALSE, clustRedDim = NULL, clustVisRes = NULL, clustVisCol = NULL) {
  observeEvent(input$closeDropDownFS1, {
    session$sendCustomMessage("close_dropDownFS", "")
  })
  ns <- session$ns
  text_input_values <- reactiveValues()
  observeEvent(inSCE$counts, {
    if (!is.null(inSCE$counts) && !is.null(colnames(colData(inSCE$counts)))) {
      options <- colnames(colData(inSCE$counts))
    } else {
      options <- NULL
    }  
    updateSelectInput(session, inputId = "ClusterLabelChoices", choices = options)
    if (is.null(inSCE$counts)) {
      print("input an SCE object")
    } else {
      count <- 0
      output$textFieldsContainer <- renderUI({
        selectInputs <- lapply(unique(colData(inSCE$counts)[[input$ClusterLabelChoices]]), function(factor) {
          local({
            div(
              style = "display: flex; align-items: center;",
              tags$label(
                style = "margin-right: 10px;",
                factor
              ),
              textInput(
                inputId = ns(paste0("textField", factor)),
                label = NULL,
                value = text_input_values[[paste0("textField", factor)]],  # Set initial value from reactiveValues
                width = "300px" 
              )
            )
          })
        })
      })
      tagList(selectInputs)
    }
    if (!is.null(inSCE$counts) && !is.null(colnames(colData(inSCE$counts)))) {
      options <- colnames(colData(inSCE$counts))
    } else {
      options <- NULL
    }  
    updateSelectInput(session, inputId = "ClusterLabelChoices", choices = options)
  })
  output$outputText <- renderText({
    "Note: If this field is left blank, the original cluster label will be updated."
  })
#  observeEvent(input$updatePlotFS1, {
#    factors <- integer()
#    text_values <- character()
#    for (i in unique(colData(inSCE$counts)[[input$ClusterLabelChoices]])) {
#      input_id <- paste0("textField", i)
#      input_value <- input[[input_id]]
#      if (!is.null(input_value) && nchar(input_value) > 0) {
#        factors <- c(factors, i)
#        text_values <- c(text_values, input_value)
#        text_input_values[[input_id]] <- input_value
#      }
#    }
#    newClusterName <- input$NewClusterLabel
#    if (newClusterName == "") {
#      inSCE$counts <- renameClusters(inSCE$counts, input$ClusterLabelChoices, factors, text_values)
#      print((colData(inSCE$counts)[,input$ClusterLabelChoices]))
#    } 
#    else {
#      inSCE$counts <- renameClusters(inSCE$counts, input$ClusterLabelChoices, factors, text_values, newClusterName)
#      updateTextInput(session, "NewClusterLabel", value = "")
#    }
#  })
  
  observeEvent(input$clustPlot1, {
    factors <- integer()
    text_values <- character()
    for (i in unique(colData(inSCE$counts)[[input$ClusterLabelChoices]])) {
      input_id <- paste0("textField", i)
      input_value <- input[[input_id]]
      if (!is.null(input_value) && nchar(input_value) > 0) {
        factors <- c(factors, i)
        text_values <- c(text_values, input_value)
        text_input_values[[input_id]] <- input_value
      }
    }
    newClusterName <- input$NewClusterLabel
    if (newClusterName == "") {
      inSCE$counts <- renameClusters(inSCE$counts, input$ClusterLabelChoices, factors, text_values)
      print((colData(inSCE$counts)[,input$ClusterLabelChoices]))
    } 
    else {
      inSCE$counts <- renameClusters(inSCE$counts, input$ClusterLabelChoices, factors, text_values, newClusterName)
      updateTextInput(session, "NewClusterLabel", value = "")
    }
    if (isCluster){
      session$sendCustomMessage("button_pressed", list(message = "Button pressed"))
    }
  })
}
