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
filterTableServer <- function(input, output, session, vals, dataframe, selectPhenotype){
  
  ns <- session$ns
  r <- reactiveValues(data = NULL, 
                      heatmapFull = NULL, 
                      heatmapPlot = NULL,
                      ridgePlot = NULL,
                      violinPlot = NULL,
                      featurePlot = NULL,
                      dotPlot = NULL
                      )
  colnames(dataframe) <- gsub("\\.", "_", colnames(dataframe))
  colnamesDF <- colnames(dataframe)
  class <- NULL
  option <- NULL
  input2 <- NULL
  for(i in seq(length(colnamesDF))){
    class[i] <- paste0("class_", colnamesDF[i])
    option[i] <- paste0("option_", colnamesDF[i])
    input2[i] <- paste0("input_", colnamesDF[i])
  }
  print(class)
  # class <- c("class1",
  #            "class2",
  #            "class3",
  #            "class4",
  #            "class5",
  #            "class6",
  #            "class7")
  
  # input1 <- c("input11",
  #             "input12",
  #             "input13",
  #             "input14",
  #             "input15",
  #             "input16",
  #             "input17")
  # 
  # input2 <- c("input21",
  #             "input22",
  #             "input23",
  #             "input24",
  #             "input25",
  #             "input26",
  #             "input27")
    
    
  # class <- c("seuratFindMarkerGeneIDDiv",
  #            "seuratFindMarkerPValDiv",
  #            "seuratFindMarkerLFCDiv",
  #            "seuratFindMarkerPct1Div",
  #            "seuratFindMarkerPct2Div",
  #            "seuratFindMarkerPValAdjDiv",
  #            "seuratFindMarkerClusterDiv"
  #            )
  # 
  # input1 <- c("seuratFindMarkerGeneIDOption",
  #                  "seuratFindMarkerPValOption",
  #                  "seuratFindMarkerLFCOption",
  #                  "seuratFindMarkerPct1Option",
  #                  "seuratFindMarkerPct2Option",
  #                  "seuratFindMarkerPValAdjOption",
  #                  "seuratFindMarkerClusterOption"
  #               )
  # input2 <- c("seuratFindMarkerGeneIDInput",
  #                 "seuratFindMarkerPValInput",
  #                 "seuratFindMarkerLFCInput",
  #                 "seuratFindMarkerPct1Input",
  #                 "seuratFindMarkerPct2Input",
  #                 "seuratFindMarkerPValAdjInput",
  #                 "seuratFindMarkerClusterInput"
  #                 )
  
  lapply(1:length(colnamesDF), function(i) {
    if(is.numeric(dataframe[,i])){
      output[[paste0("b",i)]] <- renderUI({
        hidden(div(class = class[i], wellPanel(style='border:0;',
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
      output[[paste0("b",i)]] <- renderUI({
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
                                                 choices = NULL,
                                                 label = NULL,
                                                 multiple = TRUE
                                               )
        ))
        )
      })
    }
    print(i)
    # if(i==1){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                               checkboxGroupButtons(
    #                                                                 inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                                 choices = c("=", "!="),
    #                                                                 justified = TRUE,
    #                                                                 individual = TRUE,
    #                                                                 size = "s",
    #                                                                 status = "primary"
    #                                                               ),
    #                                                               selectizeInput(
    #                                                                 inputId = ns(input2[i]),
    #                                                                 choices = NULL,
    #                                                                 label = NULL,
    #                                                                 multiple = TRUE
    #                                                               )
    #     )))
    #   })
    # }
    # if(i==2){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                             checkboxGroupButtons(
    #                                                               inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                               choices = c("<", ">", "=", "<=", ">="),
    #                                                               justified = TRUE,
    #                                                               individual = TRUE,
    #                                                               size = "s",
    #                                                               status = "primary"
    #                                                             ),
    #                                                             numericInput(
    #                                                               inputId = ns(input2[i]),
    #                                                               label = NULL,
    #                                                               step = 0.001,
    #                                                               value = 0
    #                                                             )
    #     )))
    #   })
    # }
    # if(i==3){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                            checkboxGroupButtons(
    #                                                              inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                              choices = c("<", ">", "=", "<=", ">="),
    #                                                              justified = TRUE,
    #                                                              individual = TRUE,
    #                                                              size = "s",
    #                                                              status = "primary"
    #                                                            ),
    #                                                            numericInput(
    #                                                              inputId = ns(input2[i]),
    #                                                              label = NULL,
    #                                                              step = 0.001,
    #                                                              value = 0
    #                                                            )
    #     )))
    #   })
    # }
    # if(i==4){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                             checkboxGroupButtons(
    #                                                               inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                               choices = c("<", ">", "=", "<=", ">="),
    #                                                               justified = TRUE,
    #                                                               individual = TRUE,
    #                                                               size = "s",
    #                                                               status = "primary",
    #                                                               selected = NULL
    #                                                             ),
    #                                                             numericInput(
    #                                                               inputId = ns(input2[i]),
    #                                                               label = NULL,
    #                                                               step = 0.001,
    #                                                               value = 0
    #                                                             )
    #     )))
    #   })
    # }
    # if(i==5){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                             checkboxGroupButtons(
    #                                                               inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                               choices = c("<", ">", "=", "<=", ">="),
    #                                                               justified = TRUE,
    #                                                               individual = TRUE,
    #                                                               size = "s",
    #                                                               status = "primary"
    #                                                             ),
    #                                                             numericInput(
    #                                                               inputId = ns(input2[i]),
    #                                                               label = NULL,
    #                                                               step = 0.001,
    #                                                               value = 0
    #                                                             )
    #     )))
    #   })
    # }
    # if(i==6){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                                checkboxGroupButtons(
    #                                                                  inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                                  choices = c("<", ">", "=", "<=", ">="),
    #                                                                  justified = TRUE,
    #                                                                  individual = TRUE,
    #                                                                  size = "s",
    #                                                                  status = "primary",
    #                                                                  selected = "<"
    #                                                                ),
    #                                                                numericInput(
    #                                                                  inputId = ns(input2[i]),
    #                                                                  label = NULL,
    #                                                                  step = 0.001,
    #                                                                  value = 0.05
    #                                                                )
    #     )))
    #   })
    # }
    # if(i==7){
    #   output[[paste0("b",i)]] <- renderUI({
    #     hidden(div(class = class[i], wellPanel(style='border:0;',
    #                                                                checkboxGroupButtons(
    #                                                                  inputId = ns(input1[i]), label = colnames(dataframe)[i],
    #                                                                  choices = c("=", "!="),
    #                                                                  justified = TRUE,
    #                                                                  individual = TRUE,
    #                                                                  size = "s",
    #                                                                  status = "primary"
    #                                                                ),
    #                                                                selectizeInput(
    #                                                                  inputId = ns(input2[i]),
    #                                                                  choices = NULL,
    #                                                                  label = NULL,
    #                                                                  multiple = TRUE
    #                                                                )
    #     )))
    #   })
    # }
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
                           uiOutput(ns("b1")),
                           uiOutput(ns("b2")),
                           uiOutput(ns("b3")),
                           uiOutput(ns("b4")),
                           uiOutput(ns("b5")),
                           uiOutput(ns("b6")),
                           uiOutput(ns("b7")),
                           # hidden(div(class = "seuratFindMarkerGeneIDDiv", wellPanel(style='border:0;',
                           #                                                           checkboxGroupButtons(
                           #                                                             inputId = ns("seuratFindMarkerGeneIDOption"), label = "gene.id",
                           #                                                             choices = c("=", "!="),
                           #                                                             justified = TRUE,
                           #                                                             individual = TRUE,
                           #                                                             size = "s",
                           #                                                             status = "primary"
                           #                                                           ),
                           #                                                           selectizeInput(
                           #                                                             inputId = ns("seuratFindMarkerGeneIDInput"),
                           #                                                             choices = NULL,
                           #                                                             label = NULL,
                           #                                                             multiple = TRUE
                           #                                                           )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerPValDiv", wellPanel(style='border:0;',
                           #                                                         checkboxGroupButtons(
                           #                                                           inputId = ns("seuratFindMarkerPValOption"), label = "p_val",
                           #                                                           choices = c("<", ">", "=", "<=", ">="),
                           #                                                           justified = TRUE,
                           #                                                           individual = TRUE,
                           #                                                           size = "s",
                           #                                                           status = "primary"
                           #                                                         ),
                           #                                                         numericInput(
                           #                                                           inputId = ns("seuratFindMarkerPValInput"),
                           #                                                           label = NULL,
                           #                                                           step = 0.001,
                           #                                                           value = 0
                           #                                                         )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerLFCDiv", wellPanel(style='border:0;',
                           #                                                        checkboxGroupButtons(
                           #                                                          inputId = ns("seuratFindMarkerLFCOption"), label = "avg_logFC",
                           #                                                          choices = c("<", ">", "=", "<=", ">="),
                           #                                                          justified = TRUE,
                           #                                                          individual = TRUE,
                           #                                                          size = "s",
                           #                                                          status = "primary"
                           #                                                        ),
                           #                                                        numericInput(
                           #                                                          inputId = ns("seuratFindMarkerLFCInput"),
                           #                                                          label = NULL,
                           #                                                          step = 0.001,
                           #                                                          value = 0
                           #                                                        )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerPct1Div", wellPanel(style='border:0;',
                           #                                                         checkboxGroupButtons(
                           #                                                           inputId = ns("seuratFindMarkerPct1Option"), label = "pct.1",
                           #                                                           choices = c("<", ">", "=", "<=", ">="),
                           #                                                           justified = TRUE,
                           #                                                           individual = TRUE,
                           #                                                           size = "s",
                           #                                                           status = "primary",
                           #                                                           selected = NULL
                           #                                                         ),
                           #                                                         numericInput(
                           #                                                           inputId = ns("seuratFindMarkerPct1Input"),
                           #                                                           label = NULL,
                           #                                                           step = 0.001,
                           #                                                           value = 0
                           #                                                         )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerPct2Div", wellPanel(style='border:0;',
                           #                                                         checkboxGroupButtons(
                           #                                                           inputId = ns("seuratFindMarkerPct2Option"), label = "pct.2",
                           #                                                           choices = c("<", ">", "=", "<=", ">="),
                           #                                                           justified = TRUE,
                           #                                                           individual = TRUE,
                           #                                                           size = "s",
                           #                                                           status = "primary"
                           #                                                         ),
                           #                                                         numericInput(
                           #                                                           inputId = ns("seuratFindMarkerPct2Input"),
                           #                                                           label = NULL,
                           #                                                           step = 0.001,
                           #                                                           value = 0
                           #                                                         )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerPValAdjDiv", wellPanel(style='border:0;',
                           #                                                            checkboxGroupButtons(
                           #                                                              inputId = ns("seuratFindMarkerPValAdjOption"), label = "p_val_adj",
                           #                                                              choices = c("<", ">", "=", "<=", ">="),
                           #                                                              justified = TRUE,
                           #                                                              individual = TRUE,
                           #                                                              size = "s",
                           #                                                              status = "primary",
                           #                                                              selected = "<"
                           #                                                            ),
                           #                                                            numericInput(
                           #                                                              inputId = ns("seuratFindMarkerPValAdjInput"),
                           #                                                              label = NULL,
                           #                                                              step = 0.001,
                           #                                                              value = 0.05
                           #                                                            )
                           # ))),
                           # hidden(div(class = "seuratFindMarkerClusterDiv", wellPanel(style='border:0;',
                           #                                                            checkboxGroupButtons(
                           #                                                              inputId = ns("seuratFindMarkerClusterOption"), label = "cluster",
                           #                                                              choices = c("=", "!="),
                           #                                                              justified = TRUE,
                           #                                                              individual = TRUE,
                           #                                                              size = "s",
                           #                                                              status = "primary"
                           #                                                            ),
                           #                                                            selectizeInput(
                           #                                                              inputId = ns("seuratFindMarkerClusterInput"),
                           #                                                              choices = NULL,
                           #                                                              label = NULL,
                           #                                                              multiple = TRUE
                           #                                                            )
                           # ))),
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
  
  output$findMarkerHeatmapPlotFullTopText <- renderUI({
    h6(paste("Heatmap plotted across all groups against genes with adjusted p-values <", input$seuratFindMarkerPValAdjInput))
  })
  
  output$seuratFindMarkerTable <- DT::renderDataTable({
    df <- dataframe[which(dataframe$p_val_adj < 0.05, arr.ind = TRUE),]
    seuratObject <- convertSCEToSeurat(vals$counts, scaledAssay = "seuratScaledData")
    indices <- list()
    cells <- list()
    groups <- unique(colData(vals$counts)[[selectPhenotype]])
    for(i in seq(length(groups))){
      indices[[i]] <- which(colData(vals$counts)[[selectPhenotype]] == groups[i], arr.ind = TRUE)
      cells[[i]] <- colnames(vals$counts)[indices[[i]]]
      cells[[i]] <- lapply(
        X = cells[[i]],
        FUN = function(t) gsub(
          pattern = "_",
          replacement = "-",
          x = t,
          fixed = TRUE)
      )
      Idents(seuratObject, cells = cells[[i]]) <- groups[i]
    }
    
    metadata(vals$counts)$seuratMarkersSubset <- df
    selectedGeneId <- input$seuratFindMarkerGeneIDInput
    selectedCluster <- input$seuratFindMarkerClusterInput
    updateSelectizeInput(session = session,
                         inputId = ns("seuratFindMarkerGeneIDInput"),
                         selected = selectedGeneId,
                         choices = df$gene.id)
    updateSelectizeInput(session = session,
                         inputId = ns("seuratFindMarkerClusterInput"),
                         selected = selectedCluster,
                         choices = df$cluster)
    df$p_val <- format(df$p_val, nsmall = 7)
    df$p_val_adj <- format(df$p_val_adj, nsmall = 7)
    df$pct.1 <- format(df$pct.1, nsmall = 7)
    df$pct.2 <- format(df$pct.2, nsmall = 7)
    df$avg_logFC <- format(df$avg_logFC, nsmall = 7)
    #print(head(df))
    df
  }, options = list(pageLength = 6, dom = "<'top'fl>t<'bottom'ip>", stateSave = TRUE
  ))
  
  output$seuratFindMarkerActiveFilters <- renderUI({
    panel(
      checkboxGroupInput(
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = "p_val_adj < 0.05",
        choiceValues = "p_val_adj"
      )
    )
  })
  
  observeEvent(input$seuratFindMarkerFilterRun,{
    vals$options <- list()
    for(i in seq(1:7)){
      if(i==1){
        vals$options[i] <- input$seuratFindMarkerGeneIDOption
      }
      if(i==2){
        vals$options[i] <- input$seuratFindMarkerPValOption
      }
      if(i==3){
        vals$options[i] <- input$seuratFindMarkerPValAdjOption
      }
      if(i==4){
        vals$options[i] <- input$seuratFindMarkerPct1Option
      }
      if(1==5){
        vals$options[i] <- input$seuratFindMarkerPct2Option
      }
      if(1==6){
        vals$options[i] <- input$seuratFindMarkerLFCOption
      }
      if(1==7){
        vals$options[i] <- input$seuratFindMarkerClusterOption
      }
    }
    
    # vals$seuratFindMarkerPValOption <- input$seuratFindMarkerPValOption
    # vals$seuratFindMarkerPValAdjOption <- input$seuratFindMarkerPValAdjOption
    # vals$seuratFindMarkerPct1Option <- input$seuratFindMarkerPct1Option
    # vals$seuratFindMarkerPct2Option <- input$seuratFindMarkerPct2Option
    # vals$seuratFindMarkerLFCOption <- input$seuratFindMarkerLFCOption
    # vals$seuratFindMarkerGeneIDOption <- input$seuratFindMarkerGeneIDOption
    # vals$seuratFindMarkerClusterOption <- input$seuratFindMarkerClusterOption
    updateSeuratFindMarkerTable()
  })
  
  updateSeuratFindMarkerTable <- function(){
    df <- NULL 
    parameters <- NULL
    p_val_operators <- ""
    if(!is.null(vals$options[2])){
      p_val_operators <- paste0(vals$options[2], collapse = "")
    }
    lfc_operators <- ""
    if(!is.null(vals$options[6])){
      lfc_operators <- paste0(vals$options[6], collapse = "")
    }
    pct1_operators <- ""
    if(!is.null(vals$options[4])){
      pct1_operators <- paste0(vals$options[4], collapse = "")
    }
    pct2_operators <- ""
    if(!is.null(vals$options[5])){
      pct2_operators <- paste0(vals$options[5], collapse = "")
    }
    p_val_adj_operators <- ""
    if(!is.null(vals$options[3])){
      p_val_adj_operators <- paste0(vals$options[3], collapse = "")
    }
    

    
    if(p_val_operators == ""
       && lfc_operators == ""
       && pct1_operators == ""
       && pct2_operators == ""
       && p_val_adj_operators == ""){
      df <- dataframe
      if(!is.null(vals$options[1])
         || !is.null(vals$options[7])){
        if(!is.null(vals$options[1])){
          if(vals$options[1] == "="){
            df <- df[input$seuratFindMarkerGeneIDInput,]
          }
          else if(vals$options[1] == "!="){
            df <- df[-match(input$seuratFindMarkerGeneIDInput, df$gene.id),]
          }
        }
        if(!is.null(vals$options[7])){
          if(vals$options[7] == "="){
            df <- df[which(input$seuratFindMarkerClusterInput == df$cluster),]
          }
          else if(vals$options[7] == "!="){
            df <- df[-which(input$seuratFindMarkerClusterInput == df$cluster),]
          }
        }
      }
    }
    else{
      allOperators <- c(
        "",
        p_val_operators,
        lfc_operators,
        pct1_operators,
        pct2_operators,
        p_val_adj_operators
      )
      
      allValues <- c(
        "",
        input$seuratFindMarkerPValInput,
        input$seuratFindMarkerLFCInput,
        input$seuratFindMarkerPct1Input,
        input$seuratFindMarkerPct2Input,
        input$seuratFindMarkerPValAdjInput
      )
      
      parameters <- list()
      for(i in seq(length(1:6))){
        if(allOperators[i] != ""){
          parameters$operators <- c(parameters$operators, allOperators[i])
          parameters$values <- c(parameters$values, allValues[i])
          parameters$cols <- c(parameters$cols, colnamesDF[i])
        }
      }
      parameters$operators <- na.omit(parameters$operators)
      parameters$values <- na.omit(parameters$values)
      parameters$cols <- na.omit(parameters$cols)
      
      df <- dataframe
      if(!is.null(vals$options[1])
         || !is.null(vals$options[7])){
        if(!is.null(vals$options[1])){
          if(vals$options[1] == "="){
            df <- df[input$seuratFindMarkerGeneIDInput,]
          }
          else if(vals$options[1] == "!="){
            df <- df[-match(input$seuratFindMarkerGeneIDInput, df$gene.id),]
          }
        }
        if(!is.null(vals$options[7])){
          if(vals$options[7] == "="){
            df <- df[which(input$seuratFindMarkerClusterInput == df$cluster),]
          }
          else if(vals$options[7] == "!="){
            df <- df[-which(input$seuratFindMarkerClusterInput == df$cluster),]
          }
        }
      }
      print("parameters")
      print(parameters)
      df <- singleCellTK:::.filterDF(df = df,
                                     operators = parameters$operators,
                                     cols = parameters$cols,
                                     values = parameters$values)
    }
    
    seuratObject <- convertSCEToSeurat(vals$counts, scaledAssay = "seuratScaledData")
    indices <- list()
    cells <- list()
    groups <- unique(colData(vals$counts)[[selectPhenotype]])
    for(i in seq(length(groups))){
      indices[[i]] <- which(colData(vals$counts)[[selectPhenotype]] == groups[i], arr.ind = TRUE)
      cells[[i]] <- colnames(vals$counts)[indices[[i]]]
      cells[[i]] <- lapply(
        X = cells[[i]],
        FUN = function(t) gsub(
          pattern = "_",
          replacement = "-",
          x = t,
          fixed = TRUE)
      )
      Idents(seuratObject, cells = cells[[i]]) <- groups[i]
    }
    
    r$heatmapFull <- DoHeatmap(seuratObject, features = df$gene.id)
    # output$findMarkerHeatmapPlotFull <- renderPlot({
    #   DoHeatmap(seuratObject, features = df$gene.id)
    # })
    
    output$findMarkerHeatmapPlotFullTopText <- renderUI({
      h6(paste("Heatmap plotted across all groups against genes with adjusted p-values <", input$seuratFindMarkerPValAdjInput))
    })
    
    output$seuratFindMarkerTable <- DT::renderDataTable({
      metadata(vals$counts)$seuratMarkersSubset <- df
      df$p_val <- format(df$p_val, nsmall = 7)
      df$p_val_adj <- format(df$p_val_adj, nsmall = 7)
      df$pct.1 <- format(df$pct.1, nsmall = 7)
      df$pct.2 <- format(df$pct.2, nsmall = 7)
      df$avg_logFC <- format(df$avg_logFC, nsmall = 7)
      df
    }, options = list(pageLength = 6, dom = "<'top'fl>t<'bottom'ip>", stateSave = TRUE
    ))
    
    activeFilterString <- list()
    choiceValuesFilter <- list()
    if(!is.null(parameters)){
      for(i in seq(length(parameters$cols))){
        activeFilterString[i] <- paste(parameters$cols[i], parameters$operators[i], parameters$values[i])
        choiceValuesFilter[i] <- paste(parameters$cols[i])
      }
      
      if(!is.null(vals$options[1])){
        activeFilterString <- append(activeFilterString, paste("gene.id", vals$options[1], paste(input$seuratFindMarkerGeneIDInput, collapse = ", ")))
        choiceValuesFilter <- append(choiceValuesFilter, "gene.id")
      }
      
      if(!is.null(vals$options[7])){
        activeFilterString <- append(activeFilterString, paste("cluster", vals$options[7], paste(input$seuratFindMarkerClusterInput, collapse = ", ")))
        choiceValuesFilter <- append(choiceValuesFilter, "cluster")
      }
      
      vals$activeFilterString <- as.character(activeFilterString)
      vals$choiceValuesFilter <- as.character(choiceValuesFilter)
      output$seuratFindMarkerActiveFilters <- renderUI({
        panel(
          # HTML(activeFilterString),
          checkboxGroupInput(
            inputId = ns("checkboxFiltersToRemove"),
            label = NULL,
            choiceNames = vals$activeFilterString,
            choiceValues = vals$choiceValuesFilter
          )
        )
      })
    }
    else{
      output$seuratFindMarkerActiveFilters <- renderUI({
        panel(
          HTML(paste("<span style='color:red'>No active filters!</span>")),
        )
      })
    }
    
    r$data <- df
  }
  
  seuratfindMarkerTableObserve <- observe(suspended = F, {
                                            input$seuratFindMarkerTable_rows_selected
                                            isolate({
                                              if(!is.null(input$seuratFindMarkerTable_rows_selected)){
                                                df <- metadata(vals$counts)$seuratMarkersSubset[input$seuratFindMarkerTable_rows_selected,]
                                                seuratObject <- convertSCEToSeurat(vals$counts, scaledAssay = "seuratScaledData")
                                        
                                                indices <- list()
                                                cells <- list()
                                                groups <- unique(colData(vals$counts)[[selectPhenotype]])
                                                for(i in seq(length(groups))){
                                                  indices[[i]] <- which(colData(vals$counts)[[selectPhenotype]] == groups[i], arr.ind = TRUE)
                                                  cells[[i]] <- colnames(vals$counts)[indices[[i]]]
                                                  cells[[i]] <- lapply(
                                                    X = cells[[i]], 
                                                    FUN = function(t) gsub(
                                                      pattern = "_", 
                                                      replacement = "-", 
                                                      x = t, 
                                                      fixed = TRUE)
                                                  )
                                                  Idents(seuratObject, cells = cells[[i]]) <- groups[i]
                                                }
                                                
                                                r$ridgePlot <- RidgePlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # output$findMarkerRidgePlot <- renderPlot({
                                                #   RidgePlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # })
                                                r$violinPlot <- VlnPlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # output$findMarkerViolinPlot <- renderPlot({
                                                #   VlnPlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # })
                                                r$featurePlot <- FeaturePlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # output$findMarkerFeaturePlot <- renderPlot({
                                                #   FeaturePlot(seuratObject, features = df$gene.id, ncol = 2)
                                                # })
                                                r$dotPlot <- DotPlot(seuratObject, features = df$gene.id)
                                                # output$findMarkerDotPlot <- renderPlot({
                                                #   DotPlot(seuratObject, features = df$gene.id)
                                                # })
                                                r$heatmapPlot <- DoHeatmap(seuratObject, features = df$gene.id)
                                                # output$findMarkerHeatmapPlot <- renderPlot({
                                                #   DoHeatmap(seuratObject, features = df$gene.id)
                                                # })
                                                
                                                updateTabsetPanel(session = session, inputId = ns("seuratFindMarkerPlotTabset"), selected = input$seuratFindMarkerPlotTabset)
                                                shinyjs::show(selector = ".seurat_findmarker_plots")
                                              }
                                              else {
                                                # removeTab(inputId = ns("seuratFindMarkerPlotTabset"), target = "Ridge Plot")
                                                # removeTab(inputId = ns("seuratFindMarkerPlotTabset"), target = "Violin Plot")
                                                # removeTab(inputId = ns("seuratFindMarkerPlotTabset"), target = "Feature Plot")
                                                # removeTab(inputId = ns("seuratFindMarkerPlotTabset"), target = "Dot Plot")
                                                # removeTab(inputId = ns("seuratFindMarkerPlotTabset"), target = "Heatmap Plot")
                                                # 
                                                # output$findMarkerRidgePlot <- NULL
                                                # output$findMarkerViolinPlot <- NULL
                                                # output$findMarkerFeaturePlot <- NULL
                                                # output$findMarkerDotPlot <- NULL
                                                # output$findMarkerHeatmapPlot <- NULL
                                                # 
                                                # appendTab(session = session.parent, inputId = ns("seuratFindMarkerPlotTabset"), tabPanel(title = "Ridge Plot",
                                                #                                                            panel(heading = "Ridge Plot",
                                                #                                                                  fluidRow(
                                                #                                                                    column(12, align = "center",
                                                #                                                                           panel(
                                                #                                                                             HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                #                                                                           )
                                                #                                                                    )
                                                #                                                                  )
                                                #                                                            )
                                                # )
                                                # )
                                                # appendTab(session = session.parent, inputId = ns("seuratFindMarkerPlotTabset"), tabPanel(title = "Violin Plot",
                                                #                                                            panel(heading = "Violin Plot",
                                                #                                                                  fluidRow(
                                                #                                                                    column(12, align = "center",
                                                #                                                                           panel(
                                                #                                                                             HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                #                                                                           )
                                                #                                                                    )
                                                #                                                                  )
                                                #                                                            )
                                                # )
                                                # )
                                                # appendTab(session = session.parent, inputId = ns("seuratFindMarkerPlotTabset"), tabPanel(title = "Feature Plot",
                                                #                                                            panel(heading = "Feature Plot",
                                                #                                                                  fluidRow(
                                                #                                                                    column(12, align = "center",
                                                #                                                                           panel(
                                                #                                                                             HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                #                                                                           )
                                                #                                                                    )
                                                #                                                                  )
                                                #                                                            )
                                                # )
                                                # )
                                                # appendTab(session = session.parent, inputId = ns("seuratFindMarkerPlotTabset"), tabPanel(title = "Dot Plot",
                                                #                                                            panel(heading = "Dot Plot",
                                                #                                                                  fluidRow(
                                                #                                                                    column(12, align = "center",
                                                #                                                                           panel(
                                                #                                                                             HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                #                                                                           )
                                                #                                                                    )
                                                #                                                                  )
                                                #                                                            )
                                                # )
                                                # )
                                                # appendTab(session = session.parent, inputId = ns("seuratFindMarkerPlotTabset"), tabPanel(title = "Heatmap Plot",
                                                #                                                            panel(heading = "Heatmap Plot",
                                                #                                                                  fluidRow(
                                                #                                                                    column(12, align = "center",
                                                #                                                                           panel(
                                                #                                                                             HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                #                                                                           )
                                                #                                                                    )
                                                #                                                                  )
                                                #                                                            )
                                                # )
                                                # )
                                                updateTabsetPanel(session = session, inputId = ns("seuratFindMarkerPlotTabset"), selected = input$seuratFindMarkerPlotTabset)
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
  
  observeEvent(input$seuratFindMarkerLFCOption, {
    if(length(input$seuratFindMarkerLFCOption) > 1){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerLFCOption"),
        selected = input$seuratFindMarkerLFCOption[-which(vals$options[6] == input$seuratFindMarkerLFCOption, arr.ind = TRUE)]
      )
    }
    else{
      vals$options[6] <- input$seuratFindMarkerLFCOption
    }
  })
  
  observeEvent(input$seuratFindMarkerPValOption, {
    if(length(input$seuratFindMarkerPValOption) > 1){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPValOption"),
        selected = input$seuratFindMarkerPValOption[-which(vals$options[2] == input$seuratFindMarkerPValOption, arr.ind = TRUE)]
      )
    }
    else{
      vals$options[2] <- input$seuratFindMarkerPValOption
    }
  })
  
  observeEvent(input$seuratFindMarkerPValAdjOption, {
    if(length(input$seuratFindMarkerPValAdjOption) > 1){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPValAdjOption"),
        selected = input$seuratFindMarkerPValAdjOption[-which(vals$options[3] == input$seuratFindMarkerPValAdjOption, arr.ind = TRUE)]
      )
    }
    else{
      vals$options[3] <- input$seuratFindMarkerPValAdjOption
    }
  })
  
  observeEvent(input$seuratFindMarkerPct1Option, {
    if(length(input$seuratFindMarkerPct1Option) > 1){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPct1Option"),
        selected = input$seuratFindMarkerPct1Option[-which(vals$options[4] == input$seuratFindMarkerPct1Option, arr.ind = TRUE)]
      )
    }
    else{
      vals$options[4] <- input$seuratFindMarkerPct1Option
    }
  })
  
  observeEvent(input$seuratFindMarkerPct2Option, {
    if(length(input$seuratFindMarkerPct2Option) > 1){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPct2Option"),
        selected = input$seuratFindMarkerPct2Option[-which(vals$options[5] == input$seuratFindMarkerPct2Option, arr.ind = TRUE)]
      )
    }
    else{
      vals$options[5] <- input$seuratFindMarkerPct2Option
    }
  })
  
  observeEvent(input$seuratFindMarkerRemoveAllFilters,{
    if("p_val" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPValOption"),
        selected = character(0),
        
      )
      vals$options[2] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("p_val_adj" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPValAdjOption"),
        selected = character(0)
      )
      vals$options[3] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("pct.1" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPct1Option"),
        selected = character(0)
      )
      vals$options[4] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("pct.2" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerPct2Option"),
        selected = character(0)
      )
      vals$options[5] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("avg_logFC" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerLFCOption"),
        selected = character(0)
      )
      vals$options[6] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("gene.id" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerGeneIDOption"),
        selected = character(0)
      )
      vals$options[1] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    if("cluster" %in% input$checkboxFiltersToRemove){
      updateCheckboxGroupButtons(
        session = session,
        inputId = ns("seuratFindMarkerClusterOption"),
        selected = character(0)
      )
      vals$options[7] <- NULL
      index <- which(input$checkboxFiltersToRemove == vals$choiceValuesFilter)
      updateCheckboxGroupInput(
        session = session,
        inputId = ns("checkboxFiltersToRemove"),
        label = NULL,
        choiceNames = vals$activeFilterString[-index],
        choiceValues = vals$choiceValuesFilter[-index]
      )
    }
    
    updateSeuratFindMarkerTable()
    
  })
  
  observeEvent(input$seuratFindMarkerSelectFilter,{
    shinyjs::show(selector = paste0(".class_", input$seuratFindMarkerSelectFilter))
    for(i in seq(length(class))){
      if(class[i] != paste0("class_", input$seuratFindMarkerSelectFilter)){
        shinyjs::hide(selector = paste0(".", class[i]))
      }
    }
  })
  
  
  return(r)
}