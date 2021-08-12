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
        actionButton("savePathway", "Save Scores"),
        downloadButton("downloadPathway", "Download Results")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Plot",
           
            
           radioButtons("Choice", label = NULL, c("Select Condition(s) of interest" = 'condition', 
                                                  "Select Geneset of interest" = 'geneset'), 
                                                  inline = TRUE),
                
            
           
           conditionalPanel(
             condition = sprintf("input['%s'] == 'condition'", "Choice"),
             selectInput("pathwayPlotVar",
                         "Select Condition(s) of interest for plot:", clusterChoice,
                         multiple = TRUE)
           ),
           
           conditionalPanel(
             condition = sprintf("input['%s'] == 'geneset'", "Choice"),
             uiOutput("selectGeneSets")
           ),
           
           plotOutput("pathwayPlot", height = "600px")
           
          ),
          
          
          
          tabPanel(
            
            "Results Table",
             DT::dataTableOutput("pathwayTable"),
             
              
            
            
          )
        )
      )
    )
  )
)

