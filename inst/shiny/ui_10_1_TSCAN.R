# User Interface for TSCAN  ---
shinyPanelTSCAN <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownTSCAN', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDEG', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDEGExp', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDECluster', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDEList', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownGenePlot', function(x){
                  $('html').click();
                });"),
  
  h1("TSCAN"),
  inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
  bsCollapse(
    id = "TSCANUI", 
    open = "TSCAN Date Input",
    bsCollapsePanel(
      "Calculate Pseudotime values",
      fluidRow(
        column(
          4,
          panel(
            selectInput("TSCANassayselect", "Choose an Assay:",
                        choices = c()),
            selectInput("TSCANReddim", "Select A ReducedDim:", currreddim),
            numericInput(inputId = "seed_TSCAN",
                         label = "Seed value for reproducibility of result:",
                         value = 12345,
                         step = 1),
            selectInput("clusterName", "choose cluster column", NULL, selected = NULL),
            withBusyIndicatorUI(actionButton("TSCANRun", "Run TSCAN")),
            
          )
        ),
        column(
          8,
          fluidRow(
            column(
              12,
              dropdown(
                fluidRow(
                  column(
                    12,
                    fluidRow(actionBttn(inputId = "closeDropDownTSCAN", 
                                        label = NULL, style = "simple", 
                                        color = "danger", icon = icon("times"), 
                                        size = "xs"), 
                             align = "right"),
                    
                    tags$h3("Visualization parameter"),
                    selectInput("TSCANVisRedDim", "Select A ReducedDim for visualization:", currreddim),
                    
                    actionBttn(
                      inputId = "TSCANPlot",
                      label = "Update Plot",
                      style = "bordered",
                      color = "primary",
                      size = "sm"
                    ),
                    
                  )
                ),
                
                inputId = "dropDownTSCAN",
                icon = icon("cog"),
                status = "primary",
                circle = FALSE,
                inline = TRUE
              ),
              shinyjqui::jqui_resizable(plotOutput("TSCANPlot"))

              
            )
          )
        )
      ),
      style = "primary"
    ),
    
    bsCollapsePanel(
      "Identify Expressive genes",
      fluidRow(
        column(
          4,
          panel(
            selectInput("pathIndexx", "Select Path Index:",
                           choices = "", multiple = FALSE),
            
            numericInput(inputId = "logFcThreshold_TSCAN",
                         label = "log2fcThreshold:",
                         value = 0,
                         step = 0.01),
            selectInput("discardCluster", "Select Cluster to discard:",
                           choices = "", multiple = TRUE),
            withBusyIndicatorUI(actionButton("findExpGenes", "Identify genes")),
            
          )
        ),
        column(
          8,
          fluidRow(
            column(
              12,
              tags$div(
                  class = "TSCAN_DEG_plots",
                  fluidRow(
                    tabsetPanel(
                      tabPanel("Heatmap",
                               dropdown(
                                 fluidRow(
                                   column(
                                     12,
                                     fluidRow(actionBttn(inputId = "closeDropDownDEG", 
                                                         label = NULL, style = "simple", 
                                                         color = "danger", icon = icon("times"), 
                                                         size = "xs"), 
                                              align = "right"),
                                     
                                     tags$h3("Visualization parameter"),
                                     selectInput("heatmapPathIndex", "Select Path Index:",
                                                    choices = "", multiple = FALSE),
                                     
                                     numericInput(inputId = "topGenes",
                                                  label = "Top Genes",
                                                  value = 10,
                                                  step = 1),
                                     actionBttn(
                                       inputId = "DEGPlot",
                                       label = "Update Plot",
                                       style = "bordered",
                                       color = "primary",
                                       size = "sm"
                                     ),
                                     
                                   )
                                 ),
                                 
                                 inputId = "dropDownDEG",
                                 icon = icon("cog"),
                                 status = "primary",
                                 circle = FALSE,
                                 inline = TRUE
                               ),
                               panel(
                                 shinyjqui::jqui_resizable(plotOutput(outputId = "DEGPlot"))
                               )),
                      
                      tabPanel("Expression Plot",
                               dropdown(
                                 fluidRow(
                                   column(
                                     12,
                                     fluidRow(actionBttn(inputId = "closeDropDownDEGExp", 
                                                         label = NULL, style = "simple", 
                                                         color = "danger", icon = icon("times"), 
                                                         size = "xs"), 
                                              align = "right"),
                                     
                                     tags$h3("Visualization parameter"),
                                     selectInput("expPathIndex", "Select Path Index:",
                                                    choices = "", multiple = FALSE),
                                     
                                     radioButtons("upDownRegulation", "Visualize Up/Down regulated genes:", c("increasing", "decreasing")),
                                     
                                     actionBttn(
                                       inputId = "DEGExpPlot",
                                       label = "Update Plot",
                                       style = "bordered",
                                       color = "primary",
                                       size = "sm"
                                     ),
                                     
                                   )
                                 ),
                                 
                                 inputId = "dropDownDEGExp",
                                 icon = icon("cog"),
                                 status = "primary",
                                 circle = FALSE,
                                 inline = TRUE
                               ),
                               panel(
                                 shinyjqui::jqui_resizable(plotOutput(outputId = "DEGExpPlot"))
                               )),
                      
                    )
                  )
                )
              )
          )
        )
      ),
      style = "primary"
    ),  
        
  bsCollapsePanel(
      "Identify DE genes",
      fluidRow(
        column(
          4,
          panel(
            
            selectInput("useCluster", "Select cluster of interest:",
                        choices = "", multiple = FALSE),
            numericInput(inputId = "fdrThreshold_TSCAN",
                         label = "fdrThreshold:",
                         value = 0,
                         step = 0.01),
            withBusyIndicatorUI(actionButton("findDEGenes", "Identify DE genes")),
            
          )
        ),
        column(
          8,
          fluidRow(
            column(
              12,
              tags$div(
                class = "TSCAN_DE_plots",
                fluidRow(
                  tabsetPanel(
                    tabPanel("Visualization",
                             dropdown(
                               fluidRow(
                                 column(
                                   12,
                                   fluidRow(actionBttn(inputId = "closeDropDownDECluster", 
                                                       label = NULL, style = "simple", 
                                                       color = "danger", icon = icon("times"), 
                                                       size = "xs"), 
                                            align = "right"),
                                   
                                   tags$h3("Visualization parameter"),
                                   
                                   selectInput("DEClusterRedDimNames", "Select A ReducedDim for visualization:", currreddim),
                                   selectInput("useVisCluster", "Select cluster of interest:",
                                               choices = "", multiple = FALSE),
                                  uiOutput("clusterPathIndex"),
                                   actionBttn(
                                     inputId = "DEClusterPlot",
                                     label = "Update Plot",
                                     style = "bordered",
                                     color = "primary",
                                     size = "sm"
                                   ),
                                   
                                 )
                               ),
                               
                               inputId = "dropDownDECluster",
                               icon = icon("cog"),
                               status = "primary",
                               circle = FALSE,
                               inline = TRUE
                             ),
                             panel(
                               shinyjqui::jqui_resizable(plotOutput(outputId = "DEClusterPlot"))
                             )),
                    
                    
                    tabPanel("List of genes",
                             dropdown(
                               fluidRow(
                                 column(
                                   12,
                                   fluidRow(actionBttn(inputId = "closeDropDownDEList", 
                                                       label = NULL, style = "simple", 
                                                       color = "danger", icon = icon("times"), 
                                                       size = "xs"), 
                                            align = "right"),
                                   
                                   tags$h3("List parameters"),
                                   
                                   selectInput("useListCluster", "Select cluster of interest:",
                                               choices = "", multiple = FALSE),
                                   uiOutput("clusterListPathIndex"),
                                   
                                   actionBttn(
                                     inputId = "DEClusterListPlot",
                                     label = "Generate List",
                                     style = "bordered",
                                     color = "primary",
                                     size = "sm"
                                   ),
                                   
                                 )
                               ),
                               
                               inputId = "dropDownDEList",
                               icon = icon("cog"),
                               status = "primary",
                               circle = FALSE,
                               inline = TRUE
                             ),
                             panel(
                               DT::dataTableOutput("DEClusterListPlot")
                             )),
                   
                    
                  )
                )
              )
            )
          )
          )
          ),
          style = "primary"
        
    ),
    
    bsCollapsePanel(
      "Visualize genes of interest",
      fluidRow(
        column(
          4,
          panel(
            textInput("geneName", "Enter Gene Name", "Gene Name"),
            selectInput("genesRedDimNames", "Select A ReducedDim for visualization:", currreddim),
            
            withBusyIndicatorUI(actionButton("runPlotGene", "Visualize gene on clusters")),
            
          )
        ),
        column(
          8,
          fluidRow(
            column(
              12,
              dropdown(
                fluidRow(
                  column(
                    12,
                    fluidRow(actionBttn(inputId = "closeDropDownGenePlot", 
                                        label = NULL, style = "simple", 
                                        color = "danger", icon = icon("times"), 
                                        size = "xs"), 
                             align = "right"),
                    
                    tags$h3("Visualization parameter"),
                    
                    selectInput("plotGenesRedDimNames", "Select A ReducedDim for visualization:", currreddim),
                    selectInput("useClusterForPlotGene", "Select cluster of interest:", choices = "", multiple = FALSE),
                    
                    
                    actionBttn(
                      inputId = "updatePlotGene",
                      label = "Update Plot",
                      style = "bordered",
                      color = "primary",
                      size = "sm"
                    ),
                    
                  )
                ),
                
                inputId = "dropDownGenePlot",
                icon = icon("cog"),
                status = "primary",
                circle = FALSE,
                inline = TRUE
              ),
             
              shinyjqui::jqui_resizable(plotOutput("updatePlotGene"))
              
              
            )
          )
        
      )
      ),
      style = "primary"
    )  
      
      
    ),
  nonLinearWorkflowUI(id = "nlw-Traj")
)
   