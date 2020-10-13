shinyPanelCellViewer <- fluidPage(tags$div(
  class = "container",
  h1("Cell Viewer"),
  radioGroupButtons(
    "viewertabs",
    choices = c("Scatter Plot", "Bar Plot", "Violin/Box Plot"),
    selected = NULL
  ),
  fluidRow(column(
    3,
    wellPanel(
      # Section 1 - Assay Settings
      actionButton("cv_button1", h4(strong("Select Coordinates"))),
      # open by default
      tags$div(
        id = "cv_collapse1",
        selectInput(
          inputId = "QuickAccess",
          label = NULL,
          choices = c("", approach_list, "Custom")
        ),
        #-+-+-+-+-+-X-Axis###################################
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Custom'", "QuickAccess"),
          h5(strong("X-Axis")),
          selectInput(
            "TypeSelect_Xaxis",
            h5("Data"),
            choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation")
          ),
          #Reduced Dimensions condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Xaxis"),
            selectizeInput(
              "ApproachSelect_Xaxis",
              label = h5("Approach"),
              choices = c(approach_list)
            ),
            selectInput("ColumnSelect_Xaxis", h5("Dimension"), choices = NULL)
          ),

          #Expression Assays condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Xaxis"),
            selectizeInput(
              "AdvancedMethodSelect_Xaxis",
              label = h5("Advanced Method"),
              choices = c(method_list)
            ),
            selectizeInput(
              "GeneSelect_Assays_Xaxis",
              label = h5("Feature"),
              choices = c(gene_list)
            )
          ),

          #Cell Annotation condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Xaxis"),
            selectizeInput(
              "AnnotationSelect_Xaxis",
              label = h5("Annotation"),
              choices = c(annotation_list)
            )
          ),

          #-+-+-+-+-+-Y-Axis###################################
          h5(strong("Y-Axis")),
          selectInput(
            "TypeSelect_Yaxis",
            h5("Data"),
            choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation")
          ),
          #Reduced Dimensions condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Yaxis"),
            selectizeInput(
              "ApproachSelect_Yaxis",
              label = h5("Approach"),
              choices = c(approach_list)
            ),
            selectInput("ColumnSelect_Yaxis", h5("Dimension"), choices = NULL)
          ),

          #Expression Assays condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Yaxis"),
            selectizeInput(
              "AdvancedMethodSelect_Yaxis",
              label = h5("Advanced Method"),
              choices = c(method_list)
            ),
            selectizeInput(
              "GeneSelect_Assays_Yaxis",
              label = h5("Feature"),
              choices = c(gene_list)
            )
          ),

          #Cell Annotation condition
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Yaxis"),
            selectizeInput(
              "AnnotationSelect_Yaxis",
              label = h5("Annotation"),
              choices = c(annotation_list)
            )
          )

        )
      ),

      #-+-+-+-+-+-colorby part1###################################
      tags$hr(),
      #Select Color by Data
      # Section 1 - Assay Settings
      actionButton("cv_button2", h4(strong("Color"))),

      # open by default
      tags$div(
        id = "cv_collapse2",
        selectInput(
          'TypeSelect_Colorby', h5(strong('Color By')), choices = c(
            "Single Color",
            "Reduced Dimensions",
            "Expression Assays",
            "Cell Annotation"
          ),
        ),
        # Single Color condition
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Single Color'", "TypeSelect_Colorby"),
          colourInput("Col", "", "purple", palette = 'limited')
        ),
        #Reduced Dimensions condition
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Colorby"),
          selectizeInput(
            "ApproachSelect_Colorby",
            label = h5("Approach"),
            choices = c(approach_list)
          ),
          selectInput("ColumnSelect_Colorby", h5("Dimension"), choices = NULL)
        ),

        #Expression Assays condition
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Colorby"),
          selectizeInput(
            "AdvancedMethodSelect_Colorby",
            label = h5("Advanced Method"),
            choices = c(method_list)
          ),
          selectizeInput(
            "GeneSelect_Assays_Colorby",
            label = h5("Feature"),
            choices = c(gene_list)
          )
        ),

        #Cell Annotation condition
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Colorby"),
          selectizeInput(
            "AnnotationSelect_Colorby",
            label = h5("Annotation"),
            choices = c(annotation_list)
          )
        ),

        # && output.hide_typebtns == 'show'
        #-+-+-+-+-+-colorby part2###################################
        conditionalPanel(
          condition = sprintf("input['TypeSelect_Colorby'] != 'Single Color' && output.hide_typebtns == 'show'"),
          radioButtons(
            "SelectColorType",
            label = NULL,
            choices = c("Categorical", "Continuous"),
          ),
        ),
        hidden(
          tags$div(
            id = "continuousColorConditional",
            colourInput("highColor", "High Color", "blue", "background", "limited"),
            colourInput("midColor", "Middle Color", "#666666", "background", "limited"),
            colourInput("lowColor", "Low Color", "white", "background", "limited")
          )
        ),
        hidden(tags$div(
          id = "categoricalColorConditional",
          selectizeInput(
            "colorTheme",
            label = h5("Color Scheme"),
            c("ggplot", "celda", "random")
          ),
          uiOutput("categoricalColorUI")
        )), 
        conditionalPanel(
          id="binningConditional",
          condition = sprintf("input['%s'] == 'Continuous' && input['%s'] != 'Single Color'", "SelectColorType", "TypeSelect_Colorby"),
            #sprintf("input['%s'] == 'Continuous'", "SelectColorType"),
          #checkboxInput("checkColorbinning", "Perform Binning", value = FALSE)
          h5(style="display: inline-block; margin-top: 0px; margin-bottom: 20px","Perform Binning"),
          prettyToggle(
            inputId = "checkColorbinning",
            label_on = "Yes",
            label_off = "No",
            value = FALSE,
            inline = TRUE,
            shape = "curve",
            width = "100%"
          )
          #switchInput(
          #  inputId = "checkColorbinning",
          #  onLabel = "Yes",
          #  offLabel = "No",
          #  value=FALSE,
          #  inline = TRUE,
          #  #shape = "curve"
          #  labelWidth = "50px",
          #  width = "auto"
          #)
        ),
        conditionalPanel(
          condition =  "input.checkColorbinning == 1 && output.hide_bins == 'show'",
          numericInput(
            "adjustColorbinning",
            h5("Number of Bins"),
            value = 2,
            min = 2
          )
        )
      ),
      #-+-+-+-+-+-group by###################################
      tags$hr(),
      shinyjs::useShinyjs(),
      # Section 1 - Assay Settings
      actionButton("cv_button3", h4(strong("Group By"))),
      # open by default
      tags$div(
        id = "cv_collapse3",
        selectizeInput(
          inputId = "adjustgroupby",
          label = NULL,
          choices = c("None", annotation_list)
        )
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
      conditionalPanel(
        id="violinConditional",
        condition = sprintf("input['%s'] == 'Violin/Box Plot'", "viewertabs"),
        # checkboxInput("vlnboxcheck", "Violin plot", value = FALSE),
        h5(style="display: inline-block; margin-top: 0px; margin-bottom: 20px","Use Violin Plot"),
        switchInput(
          inputId = "vlnboxcheck",
          onLabel = "Yes",
          offLabel = "No",
          value=FALSE,
          size="mini",
          inline = TRUE
        )
      ),
      actionButton("runCellViewer", "Plot")
    )
  ), #sidebarPanel_end
  #-+-+-+-+-+-mainPanel#################################
  column(
    9,
    wellPanel(
      plotlyOutput("scatter", height = "600px") %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),

      tags$br(),
      # conditionalPanel("$('#scatter').hasClass('recalculating')",
      #                  tags$div('Your plot is loading, due to large manipulation.
      #                           This message will disappear once the plot is generated.')),
      tags$hr(),
      fluidRow(
        column(6, textInput("adjusttitle", h5(strong(
          "Title:"
        )))),
        column(6, textInput("adjustlegendtitle", h5(
          strong("Legend title:")
        ))),
        column(6, textInput("adjustxlab", h5(
          strong("X-axis label:")
        ))),
        column(6, textInput("adjustylab", h5(
          strong("Y-axis label:")
        ))),
        column(3, numericInput(
          "adjustlegendtitlesize",
          h5(strong("Legend title size:")),
          min = 1,
          max = 20,
          value = 12
        )),
        column(3, numericInput(
          "adjustlegendsize",
          h5(strong("Legend size:")),
          min = 1,
          max = 20,
          value = 10
        )),
        column(3, numericInput(
          "adjustaxissize",
          h5(strong("Axis size:")),
          min = 1,
          max = 20,
          value = 10
        )),
        column(3, numericInput(
          "adjustaxislabelsize",
          h5(strong("Axis label size:")),
          min = 1,
          max = 20,
          value = 10
        )),
        column(3, numericInput(
          "adjustalpha",
          h5(strong("Opacity:")),
          min = 0,
          max = 1,
          value = 1
        )),
        column(3, numericInput(
          "adjustsize",
          h5(strong("Dot size:")),
          min = 0.1,
          max = 0.8,
          value = 0.45
        )),
        column(3, checkboxInput(
          "adjustgridlines",
          h5(strong("Add gridlines")),
          value = FALSE,
        ))
      )
    )
  ))
))

