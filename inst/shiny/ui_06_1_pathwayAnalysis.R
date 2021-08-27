shinyPanelvam <- fluidPage(
  tags$div(
    class = "container",
    h1("Pathway Analysis"),
    #h5(tags$a(href = paste0(docs.artPath, "ui_gsva.html"),
     #         "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        uiOutput("vamAssay"),
        selectInput(
                    inputId = "pathway",
                    label = "Select Method: ",
                    choices = c(
                      "VAM" = "VAM",
                      "GSVA" = "GSVA"
                    )  
        ),
        
       uiOutput("selectPathwayGeneLists"),
       
       conditionalPanel(
         condition = "input.pathway == 'VAM'",
         selectInput(
           inputId = "vamCenterParameter",
           label = "Select Center Value: ",
           choices = c(
             "TRUE" = "TRUE",
             "FALSE" = "FALSE"
           )),
         
         selectInput(
           inputId = "vamGammaParameter",
           label = "Select Gamma Value: ",
           choices = c(
             "TRUE" = "TRUE",
             "FALSE" = "FALSE"
           )),
       ),
        
       conditionalPanel(
         condition = "input.pathway == 'GSVA'",
         #radioButtons("pathwayOutPlot", "Plot Type:", c("Heatmap", "Violin"))
         
       ),
       
        withBusyIndicatorUI(actionButton("pathwayRun", "Run")),
        tags$hr(),
        h3("Save results:"),
        downloadButton("downloadPathway", "Download Results")
      ),
      mainPanel(
        

          panel(
            heading = "Visualization"
          ),
      
         
        #h2("Visualization"),
        uiOutput("selectGeneSets"),
        
        selectInput("pathwayPlotVar",
                       "Select Condition(s) of interest to group data (OPTIONAL) : ", clusterChoice,
                       multiple = TRUE),
           
        withBusyIndicatorUI(actionButton('Plot', 'Plot')),
        shinyjqui::jqui_resizable(plotOutput("depathwayPlot", height = "600px"))   
        
      )
    )
  )
)

