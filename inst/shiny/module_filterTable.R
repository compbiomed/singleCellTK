#user-interface
filterTableUI <- function(id){
  ns <- NS(id)
  fluidPage(
    fluidRow(
      panel(heading = "Marker Genes",
            uiOutput(ns("seuratFindMarkerFilter")),
            br(),
            DT::dataTableOutput(
              outputId = ns("seuratFindMarkerTable")) %>% withSpinner(
                type = 5,
                color = "#b2b2b2"
                )
      )
    )
  )
}

#server
filterTableServer <- function(input, output, session, dataframe,
                              defaultFilterColumns = NULL,
                              defaultFilterOperators = NULL,
                              defaultFilterValues = NULL){

  ns <- session$ns
  x <- session$ns('tmp')
  moduleID <- substr(x, 1, nchar(x)-4)
  rv <- reactiveValues(data = NULL,
                      selectedRows = NULL,
                      parameters = NULL
                      )
  if(length(defaultFilterValues) != length(defaultFilterColumns)
     || length(defaultFilterOperators) != length(defaultFilterColumns)){
    stop("Using a default filter requires all default filter parameters to be equal in length!")
  }
  message("Removing '.' from column names of the input dataframe as it is not supported by filters!")
  colnames(dataframe) <- gsub("\\.", "_", colnames(dataframe))
  colnamesDF <- colnames(dataframe)
  class <- NULL
  option <- NULL
  inputFirst <- NULL
  inputSecond <- NULL
  for(i in seq(length(colnamesDF))){
    class[i] <- paste0("class_", colnamesDF[i])
    option[i] <- paste0("option_", colnamesDF[i])
    inputFirst[i] <- paste0("inputFirst_", colnamesDF[i])
    inputSecond[i] <- paste0("inputSecond_", colnamesDF[i])
  }

  lapply(1:length(colnamesDF), function(i) {
    if(is.numeric(dataframe[,i])){
      if(i == 1){
        output[[paste0("filterOutput",i)]] <- renderUI({
          div(class = class[i], wellPanel(style='border:1;',
                                                 checkboxGroupButtons(
                                                   inputId = ns(option[i]),
                                                   label = colnamesDF[i],
                                                   choices = c("<", ">", "=", "<=", ">=",
                                                               "extremes", "range"),
                                                   justified = FALSE,
                                                   individual = FALSE,
                                                   size = "xs",
                                                   status = "primary"
                                                 ),
                                          conditionalPanel(
                                            condition = paste0("input['", option[i], "'] == 'extremes'"),
                                            ns = ns,
                                            h6("values greater than (or equal to):")
                                          ),
                                          conditionalPanel(
                                            condition = paste0("input['", option[i], "'] == 'range'"),
                                            ns = ns,
                                            h6("values between:")
                                          ),
                                                 numericInput(
                                                   inputId = ns(inputFirst[i]),
                                                   label = NULL,
                                                   step = 0.001,
                                                   value = 0
                                                 ),
                                          conditionalPanel(
                                            condition = paste0("input['", option[i], "'] == 'extremes'"),
                                            ns = ns,
                                            h6("and values less than (or equal to):")
                                          ),
                                          conditionalPanel(
                                            condition = paste0("input['", option[i], "'] == 'range'"),
                                            ns = ns,
                                            h6("and:")
                                          ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'extremes'
                                                                    || input['", option[i], "'] == 'range'"),
                                                   ns = ns,
                                                   numericInput(
                                                     inputId = ns(inputSecond[i]),
                                                     label = NULL,
                                                     step = 0.001,
                                                     value = 0
                                                   )
                                                 )
          ))
        })
      }
      else{
        output[[paste0("filterOutput",i)]] <- renderUI({
          hidden(div(class = class[i], wellPanel(style='border:1;',
                                                 checkboxGroupButtons(
                                                   inputId = ns(option[i]),
                                                   label = colnamesDF[i],
                                                   choices = c("<", ">", "=", "<=", ">=",
                                                               "extremes", "range"),
                                                   justified = FALSE,
                                                   individual = FALSE,
                                                   size = "xs",
                                                   status = "primary"
                                                 ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'extremes'"),
                                                   ns = ns,
                                                   h6("values greater than (or equal to):")
                                                 ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'range'"),
                                                   ns = ns,
                                                   h6("values between:")
                                                 ),
                                                 numericInput(
                                                   inputId = ns(inputFirst[i]),
                                                   label = NULL,
                                                   step = 0.001,
                                                   value = 0
                                                 ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'extremes'"),
                                                   ns = ns,
                                                   h6("and values less than (or equal to):")
                                                 ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'range'"),
                                                   ns = ns,
                                                   h6("and:")
                                                 ),
                                                 conditionalPanel(
                                                   condition = paste0("input['", option[i], "'] == 'extremes'
                                                                    || input['", option[i], "'] == 'range'"),
                                                   ns = ns,
                                                   numericInput(
                                                     inputId = ns(inputSecond[i]),
                                                     label = NULL,
                                                     step = 0.001,
                                                     value = 0
                                                   )
                                                 )
          )))
        })
      }
    }
    else if(is.character(dataframe[,i])
            || is.factor(dataframe[,i])){
      if(i == 1){
        output[[paste0("filterOutput",i)]] <- renderUI({
            div(class = class[i], wellPanel(style='border:0;',
                                            checkboxGroupButtons(
                                              inputId = ns(option[i]), label = colnamesDF[i],
                                              choices = c("=", "!="),
                                              justified = FALSE,
                                              individual = FALSE,
                                              size = "xs",
                                              status = "primary"
                                            ),
                                            selectizeInput(
                                              inputId = ns(inputFirst[i]),
                                              choices = unique(dataframe[, colnamesDF[i]]),
                                              label = NULL,
                                              multiple = TRUE
                                            )
            ))
        })
      }
      else{
        output[[paste0("filterOutput",i)]] <- renderUI({
          hidden(
            div(class = class[i], wellPanel(style='border:0;',
                                            checkboxGroupButtons(
                                              inputId = ns(option[i]), label = colnamesDF[i],
                                              choices = c("=", "!="),
                                              justified = FALSE,
                                              individual = FALSE,
                                              size = "xs",
                                              status = "primary"
                                            ),
                                            selectizeInput(
                                              inputId = ns(inputFirst[i]),
                                              choices = unique(dataframe[, colnamesDF[i]]),
                                              label = NULL,
                                              multiple = TRUE
                                            )
            ))
          )
        })
      }
    }
  })

  output$seuratFindMarkerFilter <- renderUI({
    fluidPage(
      fluidRow(
        h6("You can view the marker genes in the table below and apply custom filters to filter the table accordingly. A joint heatmap for all the marker genes available in the table is plotted underneath the table. Additional visualizations are plotted for select genes which can be selected by clicking on the rows of the table.")
      ),
      br(),
                     panel(heading = "Active Filters",
                           uiOutput(ns("seuratFindMarkerActiveFilters")),
                           br(),
                           dropdownButton(
                             fluidRow(
                               panel(
                               column(12,
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
                               )
                             ),
                             inputId = ns("addFilterDropdown"),
                             label = "Add Filter",
                             circle = FALSE,
                             inline = TRUE
                           ),
                        if(all(rv$parameters$operators == "NULL")){
                          disabled(actionButton(
                            inputId = ns("seuratFindMarkerRemoveAllFilters"),
                            label = "Remove Filter"
                          )
                          )
                        }
                        else{
                          actionButton(
                            inputId = ns("seuratFindMarkerRemoveAllFilters"),
                            label = "Remove Filter"
                          )
                        }
                     ),
    )
  })

  if(!is.null(defaultFilterColumns)){
    for(i in seq(length(colnamesDF))){
      rv$parameters$operators[i] <- "NULL"
      rv$parameters$values[i] <- "NULL"
    }
    index <- match(defaultFilterColumns, colnamesDF)
    rv$parameters$operators[index] <- defaultFilterOperators
    rv$parameters$values[index] <- defaultFilterValues

    output$seuratFindMarkerTable <- DT::renderDataTable({
      df <- .filterDF(df = dataframe,
                      operators = rv$parameters$operators,
                      cols = colnamesDF,
                      values = rv$parameters$values)
      rv$data <- df
      rv$data
    }, extensions = 'Buttons', options = list(pageLength = 6, dom = "<'top'li>t<'bottom'Bp>", stateSave = TRUE,
                                              buttons = list(
                                                list(
                                                  extend = "collection",
                                                  text = 'Export',
                                                  action = DT::JS(paste0("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('", moduleID,"-export', true, {priority: 'event'});}"))
                                                )
                                              )
    ))

    activeFilters <- list()
    activeFiltersValues <- list()
    if(!is.null(rv$parameters)){
      for(i in seq(length(colnamesDF))){
        if(rv$parameters$operators[i] != 'NULL'){
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
  }
    else{
      output$seuratFindMarkerTable <- DT::renderDataTable({
        rv$data <- dataframe
        rv$data
      }, extensions = 'Buttons', options = list(pageLength = 6, dom = "<'top'li>t<'bottom'Bp>", stateSave = TRUE,
                                                buttons = list(
                                                  list(
                                                    extend = "collection",
                                                    text = 'Export',
                                                    action = DT::JS(paste0("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('", moduleID,"-export', true, {priority: 'event'});}"))
                                                  )
                                                )
      ))

      output$seuratFindMarkerActiveFilters <- renderUI({
        panel(
          HTML(paste("<span style='color:red'>No active filters!</span>")),
        )
      })
    }

  observeEvent(input$seuratFindMarkerFilterRun,{
    #update table
    updateSeuratFindMarkerTable()

    updateSelectInput(
      session = session,
      inputId = "seuratFindMarkerSelectFilter",
      selected = colnamesDF[1]
    )
    #reset inputs
    lapply(1:length(colnamesDF), function(i) {
      if(is.numeric(dataframe[,i])){
        if(input$seuratFindMarkerSelectFilter == colnamesDF[i]){
          output[[paste0("filterOutput",i)]] <- renderUI({
            div(class = class[i], wellPanel(style='border:1;',
                                                   checkboxGroupButtons(
                                                     inputId = ns(option[i]),
                                                     label = colnamesDF[i],
                                                     choices = c("<", ">", "=", "<=", ">=",
                                                                 "extremes", "range"),
                                                     justified = FALSE,
                                                     individual = FALSE,
                                                     size = "xs",
                                                     status = "primary"
                                                   ),
                                            conditionalPanel(
                                              condition = paste0("input['", option[i], "'] == 'extremes'"),
                                              ns = ns,
                                              h6("values greater than (or equal to):")
                                            ),
                                            conditionalPanel(
                                              condition = paste0("input['", option[i], "'] == 'range'"),
                                              ns = ns,
                                              h6("values between:")
                                            ),
                                                   numericInput(
                                                     inputId = ns(inputFirst[i]),
                                                     label = NULL,
                                                     step = 0.001,
                                                     value = 0
                                                   ),
                                            conditionalPanel(
                                              condition = paste0("input['", option[i], "'] == 'extremes'"),
                                              ns = ns,
                                              h6("and values less than (or equal to):")
                                            ),
                                            conditionalPanel(
                                              condition = paste0("input['", option[i], "'] == 'range'"),
                                              ns = ns,
                                              h6("and:")
                                            ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'extremes'
                                                                    || input['", option[i], "'] == 'range'"),
                                                     ns = ns,
                                                     numericInput(
                                                       inputId = ns(inputSecond[i]),
                                                       label = NULL,
                                                       step = 0.001,
                                                       value = 0
                                                     )
                                                   )
            ))
          })
        }
        else{
          output[[paste0("filterOutput",i)]] <- renderUI({
            hidden(div(class = class[i], wellPanel(style='border:1;',
                                                   checkboxGroupButtons(
                                                     inputId = ns(option[i]),
                                                     label = colnamesDF[i],
                                                     choices = c("<", ">", "=", "<=", ">=",
                                                                 "extremes", "range"),
                                                     justified = FALSE,
                                                     individual = FALSE,
                                                     size = "xs",
                                                     status = "primary"
                                                   ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'extremes'"),
                                                     ns = ns,
                                                     h6("values greater than (or equal to):")
                                                   ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'range'"),
                                                     ns = ns,
                                                     h6("values between:")
                                                   ),
                                                   numericInput(
                                                     inputId = ns(inputFirst[i]),
                                                     label = NULL,
                                                     step = 0.001,
                                                     value = 0
                                                   ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'extremes'"),
                                                     ns = ns,
                                                     h6("and values less than (or equal to):")
                                                   ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'range'"),
                                                     ns = ns,
                                                     h6("and:")
                                                   ),
                                                   conditionalPanel(
                                                     condition = paste0("input['", option[i], "'] == 'extremes'
                                                                    || input['", option[i], "'] == 'range'"),
                                                     ns = ns,
                                                     numericInput(
                                                       inputId = ns(inputSecond[i]),
                                                       label = NULL,
                                                       step = 0.001,
                                                       value = 0
                                                     )
                                                   )
            )))
          })
        }
      }
      else if(is.character(dataframe[,i])
              || is.factor(dataframe[,i])){
        if(input$seuratFindMarkerSelectFilter == colnamesDF[i]){
          output[[paste0("filterOutput",i)]] <- renderUI({
              div(class = class[i], wellPanel(style='border:0;',
                                              checkboxGroupButtons(
                                                inputId = ns(option[i]), label = colnamesDF[i],
                                                choices = c("=", "!="),
                                                justified = FALSE,
                                                individual = FALSE,
                                                size = "xs",
                                                status = "primary"
                                              ),
                                              selectizeInput(
                                                inputId = ns(inputFirst[i]),
                                                choices = unique(dataframe[, colnamesDF[i]]),
                                                label = NULL,
                                                multiple = TRUE
                                              )
              ))
          })
        }
        else{
          output[[paste0("filterOutput",i)]] <- renderUI({
            hidden(
              div(class = class[i], wellPanel(style='border:0;',
                                              checkboxGroupButtons(
                                                inputId = ns(option[i]), label = colnamesDF[i],
                                                choices = c("=", "!="),
                                                justified = FALSE,
                                                individual = FALSE,
                                                size = "xs",
                                                status = "primary"
                                              ),
                                              selectizeInput(
                                                inputId = ns(inputFirst[i]),
                                                choices = unique(dataframe[, colnamesDF[i]]),
                                                label = NULL,
                                                multiple = TRUE
                                              )
              ))
            )
          })
        }
      }
    })

    shinyjs::enable(id = "seuratFindMarkerRemoveAllFilters")
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
      if(is.null(input[[inputFirst[i]]])){
        parameters$values[i] <- "NULL"
      }
      else{
        if(!is.null(input[[option[i]]])){
          if(input[[option[i]]] == "range" || input[[option[i]]] == "extremes"){
            parameters$values[i] <- paste0(input[[inputFirst[i]]], ",", input[[inputSecond[i]]])
          }
          else{
            parameters$values[i] <- input[[inputFirst[i]]]
          }
        }
        else{
          parameters$values[i] <- input[[inputFirst[i]]]
        }
      }
    }

    if(!is.null(rv$parameters)){
      for(i in seq(length(colnamesDF))){
        if(parameters$operators[i] != "NULL"){
          rv$parameters$operators[i] <- parameters$operators[i]
          rv$parameters$values[i] <- parameters$values[i]
        }
      }
    }
    else{
      rv$parameters <- parameters
    }

    df <- .filterDF(df = dataframe,
                                   operators = rv$parameters$operators,
                                   cols = colnamesDF,
                                   values = rv$parameters$values)





    output$seuratFindMarkerTable <- DT::renderDataTable({
      df
    }, extensions = 'Buttons', options = list(pageLength = 6, dom = "<'top'li>t<'bottom'Bp>", stateSave = TRUE,
                                              buttons = list(
                                                list(
                                                  extend = "collection",
                                                  text = 'Export',
                                                  action = DT::JS(paste0("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('", moduleID,"-export', true, {priority: 'event'});}"))
                                                )
                                              )
    ))


    activeFilters <- list()
    activeFiltersValues <- list()
    if(!is.null(rv$parameters)){
      for(i in seq(length(colnamesDF))){
        if(rv$parameters$operators[i] != 'NULL'){
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


    if(!is.null(df)){
      rv$data <- df
    }
  }

  seuratfindMarkerTableObserve <- observe(suspended = FALSE, {
                                            input$seuratFindMarkerTable_rows_selected
                                            isolate({
                                              if(!is.null(input$seuratFindMarkerTable_rows_selected)){
                                                rv$selectedRows <- input$seuratFindMarkerTable_rows_selected
                                              }
                                            })
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
    }, extensions = 'Buttons', options = list(pageLength = 6, dom = "<'top'li>t<'bottom'Bp>", stateSave = TRUE,
                                              buttons = list(
                                                list(
                                                  extend = "collection",
                                                  text = 'Export',
                                                  action = DT::JS(paste0("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('", moduleID,"-export', true, {priority: 'event'});}"))
                                                )
                                              )
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

    updateSelectInput(
      session = session,
      inputId = "seuratFindMarkerSelectFilter",
      selected = colnamesDF[index]
    )

  })

  observe({
    req(rv)
    if(all(rv$parameters$operators == "NULL")){
      shinyjs::disable(id = "seuratFindMarkerRemoveAllFilters")
      output$seuratFindMarkerActiveFilters <- renderUI({
        panel(
          HTML(paste("<span style='color:red'>No active filters!</span>")),
        )
      })
    }
    else{
      shinyjs::enable(id = "seuratFindMarkerRemoveAllFilters")
    }
  })


  observeEvent(input$export, {
    showModal(
      modalDialog(actionButton(ns("dlCSV"),"Download as CSV"),
                  br(),
                  actionButton(ns("dlPDF"),"Download as PDF"),
                  easyClose = TRUE, title = "Export Table"))
  })

  observeEvent(input$dlCSV, {
    write.csv(rv$data, file = paste0(moduleID, "-", Sys.Date(), ".csv"), row.names = TRUE)
    showNotification("Table saved in working directory as", paste0(moduleID, "-", Sys.Date(), ".csv"), duration = 10)
    removeModal()
  })

  observeEvent(input$dlPDF, {
    df <- rv$data
    dim(df)
    maxrow = 35
    npages = ceiling(nrow(df)/maxrow)
    pdf(paste0(moduleID, "-", Sys.Date(), ".pdf"), height = 11, width = 8.5)
    idx = seq(1, maxrow)
    grid.table(df[idx,],rows = NULL)
    for(i in 2:npages){
      grid.newpage();
      if(i*maxrow <= nrow(df)){
        idx = seq(1+((i-1)*maxrow), i * maxrow)
      }
      else{
        idx = seq(1+((i-1)*maxrow), nrow(df))
      }
      grid.table(df[idx, ],rows = NULL)
    }
    dev.off()
    showNotification("Table saved in working directory as", paste0(moduleID, "-", Sys.Date(), ".pdf"), duration = 10)
    removeModal()
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
      if(operators[i] == "range" || operators[i] == "extremes"){
        splitValues <- values[[i]]
        splitValues <- strsplit(splitValues, ",")
        if(operators[i] == "range"){
          filters <- c(filters, paste0("eval(call('", ">=", "', df[['", cols[i], "']],", splitValues[[1]][1], "))"))
          filters <- c(filters, paste0("eval(call('", "<=", "', df[['", cols[i], "']],", splitValues[[1]][2], "))"))
        }
        else{
          filters <- c(filters, paste0("df[['", cols[i],"']] >= ", splitValues[[1]][1]," | df[['", cols[i],"']] <= ", splitValues[[1]][2]))
        }
      }
      else{
        if(is.na(as.numeric(values[i]))){
          values[i] <- paste0("'", values[i], "'")
        }
        filters <- c(filters, paste0("eval(call('", operators[i], "', df[['", cols[i], "']],", values[i], "))"))
      }
    }
  }
  filters <- paste(filters, collapse = ",")
  parseString <- paste0("df %>% dplyr::filter(", filters, ")")
  eval(parse(text = parseString))
}
