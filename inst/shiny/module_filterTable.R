#user-interface
filterTableUI <- function(id){
  ns <- NS(id)
  fluidPage(
    fluidRow(
      tags$div(class = "seurat_findmarker_table", panel(heading = "Marker Genes",
                                                        uiOutput(ns("seuratFindMarkerFilter")),
                                                        DT::dataTableOutput(
                                                          outputId = ns("seuratFindMarkerTable")
                                                        ) %>% withSpinner(type = 5, color = "#b2b2b2")
      )
      )
    )
  )
}

#server
filterTableServer <- function(input, output, session, dataframe){
  
  ns <- session$ns
  rv <- reactiveValues(data = NULL,
                      selectedRows = NULL,
                      parameters = NULL
                      )
  message("Removing '.' from column names of the input dataframe as it is not supported by filters!")
  colnames(dataframe) <- gsub("\\.", "_", colnames(dataframe))
  colnamesDF <- colnames(dataframe)
  class <- NULL
  option <- NULL
  option2 <- NULL
  input2 <- NULL
  input3 <- NULL
  input4 <- NULL
  input5 <- NULL
  for(i in seq(length(colnamesDF))){
    class[i] <- paste0("class_", colnamesDF[i])
    option[i] <- paste0("option_", colnamesDF[i])
    option2[i] <- paste0("option2_", colnamesDF[i])
    input2[i] <- paste0("input_", colnamesDF[i])
    input3[i] <- paste0("input3_", colnamesDF[i])
    input4[i] <- paste0("input4_", colnamesDF[i])
    input5[i] <- paste0("input5_", colnamesDF[i])
  }
  
  lapply(1:length(colnamesDF), function(i) {
    if(is.numeric(dataframe[,i])){
      output[[paste0("filterOutput",i)]] <- renderUI({
        hidden(div(class = class[i], wellPanel(style='border:0;',
                                               checkboxInput(
                                                 inputId = ns(input4[i]),
                                                 label = "Double?",
                                                 value = FALSE
                                               ),
                                               conditionalPanel(
                                                 condition = paste0("input['", input4[i], "']"),
                                                 ns = ns,
                                                 radioGroupButtons(
                                                   inputId = ns(input3[i]), label = NULL,
                                                   choices = c("OR", "AND"),
                                                   justified = TRUE,
                                                   individual = TRUE,
                                                   size = "s",
                                                   status = "primary"
                                                 ),
                                                 checkboxGroupButtons(
                                                   inputId = ns(option2[i]), label = colnamesDF[i],
                                                   choices = c("<", ">", "=", "<=", ">="),
                                                   justified = TRUE,
                                                   individual = TRUE,
                                                   size = "s",
                                                   status = "primary"
                                                 ),
                                                 numericInput(
                                                   inputId = ns(input5[i]),
                                                   label = NULL,
                                                   step = 0.001,
                                                   value = 0
                                                 )
                                               ),
                                               checkboxGroupButtons(
                                                 inputId = ns(option[i]), label = colnamesDF[i],
                                                 choices = c("<", ">", "=", "<=", ">="),
                                                 justified = TRUE,
                                                 individual = TRUE,
                                                 size = "s",
                                                 status = "primary"
                                               ),
                                               numericInput(
                                                 inputId = ns(input2[i]),
                                                 label = NULL,
                                                 step = 0.001,
                                                 value = 0
                                               )
        )))
      })
    }
    else if(is.character(dataframe[,i])
            || is.factor(dataframe[,i])){
      output[[paste0("filterOutput",i)]] <- renderUI({
        hidden(
          div(class = class[i], wellPanel(style='border:0;',
                                               checkboxGroupButtons(
                                                 inputId = ns(option[i]), label = colnamesDF[i],
                                                 choices = c("=", "!="),
                                                 justified = TRUE,
                                                 individual = TRUE,
                                                 size = "s",
                                                 status = "primary"
                                               ),
                                               selectizeInput(
                                                 inputId = ns(input2[i]),
                                                 choices = unique(dataframe[, colnamesDF[i]]),
                                                 label = NULL,
                                                 multiple = TRUE
                                               )
        ))
        )
      })
    }
  })
  
  output$seuratFindMarkerFilter <- renderUI({
    fluidPage(
      fluidRow(
        h6("You can view the marker genes in the table below and apply custom filters to filter the table accordingly. A joint heatmap for all the marker genes available in the table is plotted underneath the table. Additional visualizations are plotted for select genes which can be selected by clicking on the rows of the table.")
      ),
      br(),
      fluidRow(
        column(4,offset = 0.1, style='padding:3px;', align = "center",
        ),
        column(4,offset = 0.1, style='padding:3px;', align = "center",
               radioGroupButtons(
                 inputId = ns("seuratFindMarkerFilterShowHide"), label = NULL,
                 choices = c("Show Filters" = "Show", "Hide Filters" = "Hide"),
                 justified = TRUE, status = "primary",
                 selected = "Hide",
                 size = "sm"
               )
        ),
        column(4,offset = 0.1, style='padding:3px;', align = "center",
        )
      ),
      div(class = "seuratFindMarkerShowHideDiv",
          panel(
            fluidRow(
              column(4,
                     panel(heading = "Options",
                           selectInput(
                             inputId = ns("seuratFindMarkerSelectFilter"),
                             label = "Select column to filter:",
                             choices = colnamesDF
                           ),
                           lapply(1:length(colnamesDF), function(i) {
                             uiOutput(ns(paste0("filterOutput", i)))
                           }),
                           actionButton(
                             inputId = ns("seuratFindMarkerFilterRun"),
                             label = "Apply Filter"
                           )
                     )
              ),
              column(8,
                     panel(heading = "Active Filters",
                           uiOutput(ns("seuratFindMarkerActiveFilters")),
                           br(),
                           actionButton(
                             inputId = ns("seuratFindMarkerRemoveAllFilters"),
                             label = "Remove Filter"
                           )
                     )
              )
            )
          ),
          br()
      )
    )
  })
  
    output$seuratFindMarkerActiveFilters <- renderUI({
      panel(
        HTML(paste("<span style='color:red'>No active filters!</span>")),
      )
    })
  
  
  output$seuratFindMarkerTable <- DT::renderDataTable({
    rv$data <- dataframe
    dataframe
  }, options = list(pageLength = 6, dom = "<'top'fl>t<'bottom'ip>", stateSave = TRUE
  ))
  
  observeEvent(input$seuratFindMarkerFilterRun,{
    updateSeuratFindMarkerTable()
  })
  
  updateSeuratFindMarkerTable <- function(){
    df <- NULL
    parameters <- list()
    parameters$operators <- list()
    parameters$values <- list()
    for(i in seq(length(colnamesDF))){
      if(is.null(input[[option[i]]])){
        parameters$operators[i] <- "NULL"
      }
      else{
        parameters$operators[i] <- input[[option[i]]]
      }
      if(is.null(input[[input2[i]]])){
        parameters$values[i] <- "NULL"
      }
      else{
        parameters$values[i] <- input[[input2[i]]]
      }
    }
    
    rv$parameters <- parameters
    
    df <- .filterDF(df = dataframe,
                                   operators = parameters$operators,
                                   cols = colnamesDF,
                                   values = parameters$values)
    

    

    
    output$seuratFindMarkerTable <- DT::renderDataTable({
      df
    }, options = list(pageLength = 6, dom = "<'top'fl>t<'bottom'ip>", stateSave = TRUE
    ))
    
    
    activeFilters <- list()
    activeFiltersValues <- list()
    if(!is.null(parameters)){
      for(i in seq(length(colnamesDF))){
        if(parameters$operators[i] != 'NULL'){
          activeFilters <- append(activeFilters, paste(colnamesDF[i], parameters$operators[i], parameters$values[i]))
          activeFiltersValues <- append(activeFiltersValues, colnamesDF[i])
        }
      }
    }

    output$seuratFindMarkerActiveFilters <- renderUI({
      panel(
              checkboxGroupInput(
                inputId = ns("checkboxFiltersToRemove"),
                label = NULL,
                choiceNames = as.character(activeFilters),
                choiceValues = as.character(activeFiltersValues)
              )
      )
    })
    
    
    if(!is.null(df)){
      rv$data <- df
    }
  }
  
  seuratfindMarkerTableObserve <- observe(suspended = F, {
                                            input$seuratFindMarkerTable_rows_selected
                                            isolate({
                                              if(!is.null(input$seuratFindMarkerTable_rows_selected)){
                                                rv$selectedRows <- input$seuratFindMarkerTable_rows_selected
                                              }
                                            })
                                          })
  
  observeEvent(input$seuratFindMarkerFilterShowHide,{
    if(input$seuratFindMarkerFilterShowHide == "Show"){
      shinyjs::show(selector = ".seuratFindMarkerShowHideDiv")
    }
    else{
      shinyjs::hide(selector = ".seuratFindMarkerShowHideDiv")
    }
  })
  
  observeEvent(input$seuratFindMarkerSelectFilter,{
    shinyjs::show(selector = paste0(".class_", input$seuratFindMarkerSelectFilter))
    for(i in seq(length(class))){
      if(class[i] != paste0("class_", input$seuratFindMarkerSelectFilter)){
        shinyjs::hide(selector = paste0(".", class[i]))
      }
    }
  })
  
  observeEvent(input$seuratFindMarkerRemoveAllFilters, {
    index <- match(input$checkboxFiltersToRemove, colnamesDF)
    rv$parameters$operators[index] <- "NULL"

    df <- .filterDF(df = dataframe,
                    operators = rv$parameters$operators,
                    cols = colnamesDF,
                    values = rv$parameters$values)
    
    rv$data <- df
    
    output$seuratFindMarkerTable <- DT::renderDataTable({
      df
    }, options = list(pageLength = 6, dom = "<'top'fl>t<'bottom'ip>", stateSave = TRUE
    ))
    
    activeFilters <- list()
    activeFiltersValues <- list()
    if(!is.null(rv$parameters)){
      for(i in seq(length(colnamesDF))){
        if(rv$parameters$operators[i] != 'NULL'){
          rv$parameters$operators[i]
          activeFilters <- append(activeFilters, paste(colnamesDF[i], rv$parameters$operators[i], rv$parameters$values[i]))
          activeFiltersValues <- append(activeFiltersValues, colnamesDF[i])
        }
      }
    }
    
    output$seuratFindMarkerActiveFilters <- renderUI({
      panel(
        checkboxGroupInput(
          inputId = ns("checkboxFiltersToRemove"),
          label = NULL,
          choiceNames = as.character(activeFilters),
          choiceValues = as.character(activeFiltersValues)
        )
      )
    })
    
  })
  
  
  return(rv)
}

.filterDF <- function(df, operators, cols, values){
  filters <- NULL
  for(i in seq(length(cols))){
    if(operators[i]!="NULL"){
      if(operators[i] == "="){
        operators[i] <- "=="
      }
      values[i] <- paste0("'", values[i], "'")
      filters <- c(filters, paste0("eval(call('", operators[i], "', df[['", cols[i], "']],", values[i], "))"))
    }
  }
  filters <- paste(filters, collapse = ",")
  parseString <- paste0("df %>% filter(", filters, ")")
  eval(parse(text = parseString))
}