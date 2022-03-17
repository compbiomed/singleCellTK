

shinyPanelvam <- fluidPage(

  tags$script("Shiny.addCustomMessageHandler('close_dropDownPathway', function(x){
                  $('html').click();
                });"),

  tags$div(
    class = "container",
    h1("Pathway Analysis"),
    h5(tags$a(href = paste0(docs.artPath, "pathwayAnalysis.html"),
              "(help)", target = "_blank")),

    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          inputId = "vamAssay", 
          label = "Select input matrix:", 
          choices = NULL, 
          selected = NULL, 
          multiple = FALSE,
          options = NULL),
        #uiOutput("vamAssay"),

        selectInput(
                    inputId = "pathway",
                    label = "Select Method: ",
                    choices = c(
                      "VAM" = "VAM",
                      "GSVA" = "GSVA"
                    )
        ),

       #uiOutput("selectPathwayGeneLists"),
       selectizeInput("PathwayGeneLists", "Select Geneset Collection(s):",
                      choices = "Import geneset before using", multiple = FALSE),

       conditionalPanel(
         condition = "input.pathway == 'VAM'",
         checkboxInput("vamCenterParameter", "Select center value",
                       value = TRUE),
         checkboxInput("vamGammaParameter", "Select gamma value", 
                       value = TRUE),
       ),

       withBusyIndicatorUI(actionButton("pathwayRun", "Run")),
       tags$hr(),
       #h3("Save results:"),
       #downloadButton("downloadPathway", "Download Results")

       ),

      mainPanel(

        dropdown(
          fluidRow(
            column(
              12,
              fluidRow(actionBttn(inputId = "closeDropDownPathway", 
                                  label = NULL, style = "simple", 
                                  color = "danger", icon = icon("times"), 
                                  size = "xs"), 
                       align = "right"),
              
              tags$h3("List of Input"),
              selectizeInput("pathwayRedDimNames", 
                             "Select Score matrix which you want to plot:", 
                             choices = NULL, 
                             multiple = FALSE),
              selectizeInput("pathwayPlotGS", "Select Geneset:",
                             choices = "", multiple = FALSE),
              selectInput(inputId = 'pathwayPlotVar',
                          label = 'Select Condition(s) of interest to group data (OPTIONAL) : ',
                          choices = clusterChoice,
                          multiple = FALSE),
              checkboxInput("pathwayPlotBoxplot", "Add box plot", 
                            value = FALSE),
              checkboxInput("pathwayPlotViolinplot", "Add violin plot", 
                            value = TRUE),
              checkboxInput("pathwayPlotDots", "Add dots", 
                            value = TRUE),
              pickerInput(
                inputId = "pathwayPlotSummary",
                label = "Select Summary parameter: ",
                choices = c("mean", "median"), 
                selected = "median"),
              
              #withBusyIndicatorUI(actionButton("Plot", "Update Plot")),
              actionBttn(
                inputId = "pathwayPlot",
                label = "Update Plot",
                style = "bordered",
                color = "primary",
                size = "sm"
              ),
              
            )
          ),
          
          #tags$style(".btn-custom {background-color: #e5e5e5; ;}"),
          #label = "Plot Options",
          inputId = "dropDownPathway",
          icon = icon("cog"),
          status = "primary",
          circle = FALSE,
          inline = TRUE
        ),
        shinyjqui::jqui_resizable(plotOutput("pathwayPlot"))
      )
    )
  )
)

