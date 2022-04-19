shinyPanelLabelCellType <- fluidPage(
  tags$div(
    class = "container",
    h1("Label Cell Type"),
    p("Currently only 'SingleR' method supported. ", style = "color: grey;"),
    h5(tags$a(href = paste0(docs.artPath, "cell_type_labeling.html"),
              "(help)", target = "_blank")),
    panel(
      fluidRow(
        column(
          width = 6,
          selectInput("ctLabelMethod", "Method", c("SingleR"), "SingleR"),
        ),
        column(
          width = 6,
          selectizeInput(
            inputId = "ctLabelAssay", 
            label = "Select input matrix:", 
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = NULL)
          #uiOutput("ctLabelAssay"),
        )
      ),
      panel(
        useShinyjs(),
        # SingleR specific panel, in case adding other methods in the future.
        p("SingleR works with a reference dataset where the cell type ",
          "labeling is given. Given a reference dataset of samples ",
          "(single-cell or bulk) with known labels, it assigns those labels ",
          "to new cells from a test dataset based on similarities in their ",
          "expression profiles.",
          style = "color: grey;"),
        fluidRow(
          column(
            width = 6,
            selectInput("ctLabelRef", "Choose a reference",
                        choices = c("Human Primary Cell Atlas Data" = "hpca",
                                    "Blueprint Encode Data" = "bpe",
                                    "Muraro Pancreas Data" = "mp",
                                    "Database Immune Cell Expression Data" = "dice",
                                    "ImmGen Data" = "immgen",
                                    "Mouse RNAseq Data" = "mouse",
                                    "Zeisel Brain Data" = "zeisel")),
            radioButtons("ctLabelFeatureType", "Feature type:",
                         c("symbol", "ensembl"), inline = TRUE),
          ),
          column(
            width = 6,
            uiOutput("ctLabelLevelUI"),
            # selectInput("ctLabelLevel", "Labeling level:",
            #             c("main", "fine", "ont"), "main"),
            radioButtons("ctLabelBy", "Label by",
                         c("Each cell", "Clusters"), "Each cell", inline = TRUE),
          )
        ),
        conditionalPanel(
          "input.ctLabelBy == 'Clusters'",
          selectInput("ctLabelByCluster", "Use cluster label:", clusterChoice)
          ),
        conditionalPanel(
          "input.ctLabelRef == 'hpca'",
          tags$a(href = "https://rdrr.io/github/LTLA/celldex/man/HumanPrimaryCellAtlasData.html",
                "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'bpe'",
          tags$a(href = "https://rdrr.io/github/LTLA/celldex/man/BlueprintEncodeData.html",
                "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'mp'",
          tags$a(href = "https://rdrr.io/github/LTLA/scRNAseq/man/MuraroPancreasData.html",
                 "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'dice'",
          tags$a(href = "https://rdrr.io/github/LTLA/celldex/man/DatabaseImmuneCellExpressionData.html",
                 "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'immgen'",
          tags$a(href = "https://rdrr.io/github/LTLA/celldex/man/ImmGenData.html",
                 "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'mouse'",
          tags$a(href = "https://rdrr.io/github/LTLA/celldex/man/MouseRNAseqData.html",
                 "Detail about the reference dataset...", target = "_blank")
        ),
        conditionalPanel(
          "input.ctLabelRef == 'zeisel'",
          tags$a(href = "https://rdrr.io/github/LTLA/scRNAseq/man/ZeiselBrainData.html",
                 "Detail about the reference dataset...", target = "_blank")
        ),
        p("If users want to work with their customized reference dataset, ",
          "please refer to SCTK's console function: runSingleR()",
          style = "color: grey;")
      ),
      withBusyIndicatorUI(actionButton("ctLabelRun", "Label"))
    )
  )
)
