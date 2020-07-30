shinyPanelGeneSets <- fluidPage(
  tags$div(
    class = "container",
    style = "margin-bottom: 10px",
    h1("Import Gene Sets"),
    tags$hr(),
    
    radioButtons("geneSetSourceChoice", label = NULL, c("Upload a GMT file" = 'gsGMTUpload',
                                                 "Select from a database" = "gsDBUpload",
                                                 "Paste in your gene set" = "gsPasteUpload")
    ),
    
    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsGMTUpload'", "geneSetSourceChoice"),
      h3("Upload a GMT file:"),
      textInput('gsCollectionNameGMT', label='Collection Name'),
      fileInput('geneSetGMT', 'Choose GMT File', accept = ".gmt")
      # shinyFilesButton('geneSetGMT', 'Choose File', 'Select a .gmt file', FALSE)
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsDBUpload'", "geneSetSourceChoice"),
      h3("Select from a database:"),
      tags$style(HTML("#geneSetDB {width:100%}")),
      checkboxGroupInput('geneSetDB', 'Check the gene sets you want to import', 
                         choices = c()),
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsPasteUpload'", "geneSetSourceChoice"),
      h3("Paste in your gene set:"),
      textInput('gsCollectionNameText', label='Create a new collection (enter collection name)'),
      shinyjs::hidden(
        tags$div(id = "gsAddToExisting",
                 h4("-OR-"),
                 selectInput("gsExisting", "Add to an existing collection", c("None")),
                 textAreaInput('geneSetText', 'Please enter values separated by new lines', width = "300px")
        )
      ),
      # shinyjs::hidden(
      #   tags$div(id = "gsAddToExisting",
      #            h4("-OR-"),
                 # selectInput("gsExisting", "Add to an existing collection", c("None")),
                 # textAreaInput('geneSetText', 'Please enter values separated by new lines', width = "300px")
      #   ),
      # )
    ),
    
    withBusyIndicatorUI(
      actionButton("uploadGS", "Upload")
    ),
    shinyjs::hidden(
      tags$div(id = "gsUploadError",
               tags$b("Please fill out all the required fields", style = "color: red;"),
      )
    ),
  )
)