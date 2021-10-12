shinyPanelGeneSets <- fluidPage(
  tags$div(
    class = "container",
    style = "margin-bottom: 10px",
    h1("Import Gene Sets"),
    h5(tags$a(href = paste0(docs.artPath, "import_genesets.html"),
              "(help)", target = "_blank")),
    tags$hr(),
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")
    ),
    wellPanel(
      h4("Existing Gene Sets:"),
      fluidRow(
        column(3, tags$b("Collection Name")),
        column(9, tags$b("Source")),
      ),
      tags$div(id = "newGSImport"),
      tags$br(),
      tags$br(),
    ),
    radioButtons("geneSetSourceChoice", label = NULL, c("Upload a GMT file" = 'gsGMTUpload',
                                                 "Select from a database" = "gsDBUpload",
                                                 "Import mitochondrial gene set" = "gsMito",
                                                 "Paste in your gene set" = "gsPasteUpload"
                                                 )
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsGMTUpload'", "geneSetSourceChoice"),
      h3("Upload a GMT file:"),
      fileInput('geneSetGMT', 'Choose GMT File', accept = ".gmt"),
      textInput('gsCollectionNameGMT', label='Collection Name'),
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsDBUpload'", "geneSetSourceChoice"),
      h3("Select from a database:"),
      tags$style(HTML("#geneSetDB {width:100%}")),
      checkboxGroupInput('geneSetDB', 'Check the gene sets you want to import',
                         choices = c()),
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsMito'", "geneSetSourceChoice"),
      h3("Import mitochondrial gene set"),
      radioButtons("geneSetMitoSpecies", "Species",
                   choices = c("human", "mouse"), selected = "human",
                   inline = TRUE),
      selectInput("geneSetMitoID", "ID Type",
                  choices = c("symbol", "entrez", "ensembl", "ensemblTranscriptID"),
                  selected = "symbol"),
      textInput("geneSetMitoName", "Collection Name", "Mito",
                placeholder = "Required"),
    ),


    conditionalPanel(
      condition = sprintf("input['%s'] == 'gsPasteUpload'", "geneSetSourceChoice"),
      h3("Paste in your gene set:"),
      textInput('gsCollectionNameText', label='Create a new collection (enter collection name)'),
      shinyjs::hidden(
        tags$div(id = "gsAddToExisting",
                 h4("-OR-"),
                 selectInput("gsExisting", "Add to an existing collection", c("None")),
        )
      ),
      textAreaInput('geneSetText', 'Please enter values separated by new lines', width = "300px")
    ),

    selectInput("gsByParam", "Location within SCE object where the gene identifiers in should be mapped.", list()),

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

