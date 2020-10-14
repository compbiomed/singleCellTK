shinyPanelCellViewer <- fluidPage(
  tags$div(
    class = "container",
           h1("Cell Viewer"),
           fluidRow(
           radioGroupButtons("viewertabs", choices = c("reducedDims Plot", "Bar Plot","Violin/Box Plot"), selected = NULL),
                    wellPanel(style = "background: floralwhite",
                      conditionalPanel(condition = sprintf("input['%s'] == 'Violin/Box Plot'", "viewertabs"),
                        checkboxInput("vlnboxcheck", "Violin plot", value = FALSE)),
                            # Section 1 - Assay Settings
                            actionButton("cv_button1", h4(strong("Select Coordinates")),style = "background: floralwhite"),
                            # open by default
                            tags$div(id = "cv_collapse1",
                            selectInput(inputId = "QuickAccess", label = NULL, choices = c("",approach_list,"Custom")),
                            #-+-+-+-+-+-X-Axis###################################
                            conditionalPanel(condition = sprintf("input['%s'] == 'Custom'", "QuickAccess"),
                            h5(strong("Select X-Axis:")),
                            selectInput("TypeSelect_Xaxis", h5("Type Of Data:"), choices = c("Reduced Dimensions","Expression Assays","Cell Annotation")),
                            #Reduced Dimensions condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Xaxis"),
                                             selectizeInput("ApproachSelect_Xaxis", label = h5("Select Approach:"), choices = c(approach_list)),
                                             selectInput("ColumnSelect_Xaxis", h5("Select Dimension:"),choices = NULL)),

                            #Expression Assays condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Xaxis"),
                                             selectizeInput("AdvancedMethodSelect_Xaxis", label = h5("Select Advanced Method:"), choices = c(method_list)),
                                             selectizeInput("GeneSelect_Assays_Xaxis", label = h5("Select Feature:"), choices = c(gene_list))),

                            #Cell Annotation condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Xaxis"),
                                             selectizeInput("AnnotationSelect_Xaxis", label = h5("Select Annotation:"), choices = c(annotation_list))),

                            #-+-+-+-+-+-Y-Axis###################################
                            tags$hr(),
                            h5(strong("Select Y-Axis:")),
                            selectInput("TypeSelect_Yaxis", h5("Type Of Data:"), choices = c("Reduced Dimensions","Expression Assays","Cell Annotation")),
                            #Reduced Dimensions condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Yaxis"),
                                             selectizeInput("ApproachSelect_Yaxis", label = h5("Select Approach:"), choices = c(approach_list)),
                                             selectInput("ColumnSelect_Yaxis", h5("Select Dimension:"),choices = NULL)),

                            #Expression Assays condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Yaxis"),
                                             selectizeInput("AdvancedMethodSelect_Yaxis", label = h5("Select Advanced Method:"), choices = c(method_list)),
                                             selectizeInput("GeneSelect_Assays_Yaxis", label = h5("Select Feature:"), choices = c(gene_list))),

                            #Cell Annotation condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Yaxis"),
                                             selectizeInput("AnnotationSelect_Yaxis", label = h5("Select Annotation:"), choices = c(annotation_list)))
                              ) #conditionPanel_end
                              ), #div_end

                            #-+-+-+-+-+-colorby part1###################################
                            tags$hr(),
                            #Select Color by Data
                            # Section 1 - Assay Settings
                            actionButton("cv_button2", h4(strong("Color")),style = "background: floralwhite"),
                            # open by default
                            tags$div(id = "cv_collapse2",
                            radioGroupButtons(inputId = "TypeSelect_Colorby", label = h5(strong("Type of Data:")), choices = c("Pick a Color","Reduced Dimensions", "Expression Assays","Cell Annotation"), direction = "horizontal"),
                            #Reduced Dimensions condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Colorby"),
                                            selectizeInput("ApproachSelect_Colorby", label = h5("-> Approach:"),
                                                           choices = c(approach_list)),
                                            selectInput("ColumnSelect_Colorby", h5("-> Dimension:"),choices = NULL)),

                            #Expression Assays condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Colorby"),
                                            selectizeInput("AdvancedMethodSelect_Colorby", label = h5("-> Advanced Method:"),
                                                           choices = c(method_list)),
                                            selectizeInput("GeneSelect_Assays_Colorby", label = h5("-> Feature:"),
                                                           choices = c(gene_list))),

                            #Cell Annotation condition
                            conditionalPanel(condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Colorby"),
                                            selectizeInput("AnnotationSelect_Colorby", label = h5("-> Annotation:"),
                                                           choices = c(annotation_list))),

                            #-+-+-+-+-+-colorby part2###################################
                            conditionalPanel(condition = sprintf("input['%s'] == 'Pick a Color'", "TypeSelect_Colorby"),
                                            colourInput("Col", h5(strong("Color Picker:")), "purple",palette = 'limited')),

                            conditionalPanel(condition = sprintf("input['%s'] != 'Pick a Color'", "TypeSelect_Colorby"),
                                            radioButtons("SelectColorType",label = NULL,choices = c("Categorical", "Continuous")),
                                            tags$hr(),

                                            conditionalPanel(condition = sprintf("input['%s'] == 'Continuous'", "SelectColorType"),
                                                             checkboxInput("checkColorbinning",h5("Perform Binning:"), value = FALSE)),

                                            conditionalPanel(condition =  "input.checkColorbinning == 1",
                                                             numericInput("adjustColorbinning", h5("Number of Bins:"), value = 2, min =2))
                                            #,


                                            #selectizeInput("adjustbrewer", h5(strong("Color Palettes:")), choices = NULL)
                            )
),
                              #-+-+-+-+-+-group by###################################
                              tags$hr(),
                              shinyjs::useShinyjs(),
                              # Section 1 - Assay Settings
                              actionButton("cv_button3", h4(strong("Group by")),style = "background: floralwhite"),
                              # open by default
                              tags$div(id = "cv_collapse3",
                                       selectizeInput(inputId = "adjustgroupby", label = NULL, choices = c("None", annotation_list))
                                #,
                                #       conditionalPanel(condition = sprintf("input['%s'] != 'None'", "adjustgroupby"),
                                #                       radioButtons("SelectValueType",label = NULL,choices = c("Categorical", "Continuous")),
                                #                       conditionalPanel(condition = sprintf("input['%s'] == 'Continuous'", "SelectValueType"),
                                #                                         checkboxInput("checkbinning",h5("Perform Binning:"), value = FALSE)),
                                #                       conditionalPanel(condition = "input.checkbinning == 1",
                                #                                         numericInput("adjustbinning", h5("Number of Bins:"),value = 2, min =2))
                                #       )
                              ),
                              tags$hr(),
                              actionButton("runCellViewer", "Plot"),
                 fluidRow(column(6,textInput("adjusttitle", h5(strong("Title:")))),
                   column(6,textInput("adjustlegendtitle", h5(strong("Legend title:")))),
                   column(6,sliderInput("adjustalpha", h5(strong("Opacity:")), min = 0, max = 1, value = 1)),
                   column(6,sliderInput("adjustsize", h5(strong("Size:")), min = 0.1, max = 0.8, value = 0.45)),
                   column(6,textInput("adjustxlab", h5(strong("X-axis label:")))),
                   column(6,textInput("adjustylab", h5(strong("Y-axis label:"))))
               )
                    )
), #fluidrow_end
      fluidRow(
             #-+-+-+-+-+-mainPanel#################################
             wellPanel(style = "background: floralwhite",
                                plotlyOutput("scatter", height = "auto") %>% withSpinner(size = 3, color="#0dc5c1", type = 8)
               )
                                # conditionalPanel("$('#scatter').hasClass('recalculating')",
                                #                  tags$div('Your plot is loading, due to large manipulation.
                                #                           This message will disappear once the plot is generated.')),
        )
)#tag_end
)#page_end

