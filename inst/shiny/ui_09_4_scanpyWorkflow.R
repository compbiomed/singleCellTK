# QC & Filtering

# Normalize

# HVG

# PCA

# Clustering

# Find Marker

# User Interface for Seurat Workflow ---
shinyPanelScanpy <- fluidPage(
    tags$script("Shiny.addCustomMessageHandler('close_dropDownScanpy', function(x){
                  $('html').click();
                });"), # close_dropDownSeuratHM
    h1("Scanpy"),
    h5(tags$a(href = paste0(docs.artPath, "cnsl_scanpy_curated_workflow.html"), #cnsl_seurat_curated_workflow.html
              "(help)", target = "_blank")),
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    conditionalPanel(
        condition = "false",
        selectInput(
            "activePanelSelectScanpy", #activePanelSelectSeurat
            label = "Active Panel:",
            choices = c("",
                        "Normalize Data",
                        "Scale Data",
                        "Highly Variable Genes",
                        "Dimensionality Reduction",
                        "tSNE/UMAP",
                        "Clustering",
                        "Find Markers",
                        "Heatmap Plot"),
            selected = ""
        )
    ),
    bsCollapse(id = "ScanpyUI", open = "QC & Filtering", # SeuratUI
            bsCollapsePanel("QC & Filtering",
                fluidRow(
                        ),
                            style = "primary"
                            ),
            bsCollapsePanel("Normalize Data",
                            fluidRow(
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Highly Variable Genes",
                            fluidRow(
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Dimensionality Reduction",
                            fluidRow(
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Clustering",
                            fluidRow(
                            ),
                            style = "primary"
            )
       ),
    nonLinearWorkflowUI(id = "nlw-scanpy") #nlw-seurat
    )