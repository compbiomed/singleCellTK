# User Interface for TSCAN  ---
shinyPanelTSCAN <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownTSCAN', function(x){
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
  
  h1("TSCAN"),
  inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
  bsCollapse(
    id = "TSCANUI", 
    open = "TSCAN Date Input",
    bsCollapsePanel(
      "Calculate Pseudotime Values",
      fluidRow(
        column(
          4,
          panel(
            
            selectInput("TSCANReddim", "reducedDim Name:", currreddim),
            numericInput(inputId = "seed_TSCAN",
                         label = "Seed value for reproducibility of result:",
                         value = 12345,
                         step = 1),
            selectInput("clusterName", "Name of Clustering Result: ", "Auto generate clusters", selected = NULL),
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
                    
                    tags$h3("Visualization Parameter"),
                    selectInput("TSCANVisRedDim", "reducedDim Name:", currreddim),
                    
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
      "Identify Genes Differentially Expressed For Path",
      fluidRow(
        column(
          4,
          panel(
            selectInput("TSCANassayselect", "Choose an assay:",
                        choices = c()),
            
            pickerInput("pathIndexx", "Select path terminal node:",
                           choices = "", multiple = FALSE),
            
            numericInput(inputId = "logFcThreshold_TSCAN",
                         label = "Log2FC greater than:",
                         value = 0,
                         step = 0.01),
            uiOutput("discardCluster"),
            #selectInput("discardCluster", "Select cluster to discard (OPTIONAL):",
             #              choices = "", multiple = TRUE),
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
                    dropdown(
                      fluidRow(
                        column(
                          12,
                          fluidRow(actionBttn(inputId = "closeDropDownDEGExp", 
                                              label = NULL, style = "simple", 
                                              color = "danger", icon = icon("times"), 
                                              size = "xs"), 
                                   align = "right"),
                          
                          tags$h3("Visualization Parameter"),
                          
                          selectInput("expPathIndex", "Select path terminal node:",
                                      choices = "", multiple = FALSE),
                          
                          numericInput(inputId = "topGenes",
                                       label = "Top Genes",
                                       value = 10,
                                       step = 1),

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
                    tabsetPanel(
                      
                      tabPanel("Heatmap",
                              
                               panel(
                                 shinyjqui::jqui_resizable(plotOutput(outputId = "heatmapPlot"))
                               )),
                      
                      tabPanel("Top Upgregulated Genes",
                              
                               panel(
                                 #shinyjqui::jqui_resizable(plotOutput(outputId = "DEGExpPlot"))
                                 shinyjqui::jqui_resizable(plotOutput(outputId = "UpregGenesPlot"))
                               )),
                      
                      tabPanel("Top Downregulated Genes",
                              
                               panel(
                                 shinyjqui::jqui_resizable(plotOutput(outputId = "DownregGenesPlot"))
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
      "Identify Genes Differentially Expressed For Branched Cluster",
      fluidRow(
        column(
          4,
          panel(
            
            selectInput("useCluster", "Select branched cluster of interest:",
                        choices = "", multiple = FALSE),
            numericInput(inputId = "fdrThreshold_TSCAN",
                         label = "FDR less than:",
                         value = 0.05,
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
                                   
                                   tags$h3("Visualization Parameter"),
                                   
                                   selectInput("DEClusterRedDimNames", "reducedDim Name:", currreddim),
                                   selectInput("useVisCluster", "Select branched cluster of interest:",
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
                                   
                                   tags$h3("List Parameters"),
                                   
                                   selectInput("useListCluster", "Select branched cluster of interest:",
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
      "Plot expression of individual genes",
      fluidRow(
        column(
          4,
          panel(
            textInput("geneName", "Enter gene name", "Gene Name"),
            selectInput("genesRedDimNames", "reducedDim Name:", currreddim),
            pickerInput("useClusterForPlotGene", "Select cluster of interest:", choices = "", multiple = TRUE, options = list(
              #`actions-box` = TRUE,
              `none-selected-text` = "All clusters selected"
            )),
            
            withBusyIndicatorUI(actionButton("runPlotGene", "Plot")),
            
          )
        ),
        column(
          8,
          fluidRow(
            column(
              12,
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
   