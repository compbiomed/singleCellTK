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
        
       withBusyIndicatorUI(actionButton("pathwayRun", "Run")),
       tags$hr(),
       #h3("Save results:"),
       #downloadButton("downloadPathway", "Download Results")
      
       ),
      mainPanel(
        
        dropdown(
          
          tags$h3("List of Input"),
          
          uiOutput("selectReduceDim"),
          
          uiOutput("selectGeneSets"),
          
          selectInput(inputId = 'pathwayPlotVar',
                      label = 'Select Condition(s) of interest to group data (OPTIONAL) : ', 
                      choices = clusterChoice,
                      multiple = FALSE),
          
          radioButtons("boxplot", "Box Plot:", 
                       c("TRUE" = "TRUE", 
                         "FALSE" = "FALSE"), inline = TRUE),
          
          radioButtons("violinplot", "Violin Plot:", 
                       c("TRUE" = "TRUE", 
                         "FALSE" = "FALSE"), inline = TRUE),
          
          pickerInput(
            inputId = "summary",
            label = "Select Summary parameter: ",
            choices = c(
              "Mean" = "mean",
              "Median" = "median"
            )),
          
          withBusyIndicatorUI(actionButton("Plot", "Update Plot")),  
          
          tags$style(".btn-custom {background-color: #e5e5e5; ;}"),
          label = "Plot Options",
          status = "custom", width = "300px", 
          animate = animateOptions(
            enter = animations$fading_entrances$fadeInLeftBig,
            exit = animations$fading_exits$fadeOutRightBig
          )
        ),
        
        
        plotOutput("pathwayPlot", height = "600px") 
          
        
      )
    )
  )
)

