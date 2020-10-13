orgpkgs <- c("Anopheles" = "org.Ag.eg.db", "Arabidopsis" = "org.At.tair.db",
             "Bovine" = "org.Bt.eg.db", "Worm" = "org.Ce.eg.db",
             "Canine" = "org.Cf.eg.db", "Fly" = "org.Dm.eg.db",
             "Zebrafish" = "org.Dr.eg.db",
             "E coli strain K12" = "org.EcK12.eg.db",
             "E coli strain Sakai" = "org.EcSakai.eg.db",
             "Chicken" = "org.Gg.eg.db", "Human" = "org.Hs.eg.db",
             "Mouse" = "org.Mm.eg.db", "Rhesus" = "org.Mmu.eg.db",
             "Malaria" = "org.Pf.plasmo.db", "Chimp" = "org.Pt.eg.db",
             "Rat" = "org.Rn.eg.db", "Yeast" = "org.Sc.sgd.db",
             "Pig" = "org.Ss.eg.db", "Xenopus" = "org.Xl.eg.db")

shinyPanelFilter <- fluidPage(
  useShinyalert(),
  wellPanel(
    h4("Select Column Filtering Criteria:"),
    fluidRow(
      column(4, tags$b("Column Name")),
      column(4, tags$b("Filter Condition")),
      column(4, tags$b("Remove"))
    ),
    tags$div(id = "newFilteringParams"),
    tags$br(),
    tags$br(),
    actionButton("addFilteringParam", "Add a Filter"),
    actionButton("clearAllFilters", "Clear Filters"),
  ),
  
  # *** Uncomment when row filtering is finalized ***
  wellPanel(
    h4("Select Row Filtering Criteria:"),
    fluidRow(
      column(4, tags$b("Assay Name")),
      column(4, tags$b("Filter Condition")),
      column(4, tags$b("Remove"))
    ),
    tags$div(id = "newRowFilteringParams"),
    tags$br(),
    tags$br(),
    actionButton("addRowFilteringParam", "Add a Filter"),
    actionButton("clearAllRowFilters", "Clear Filters"),
  ),
  withBusyIndicatorUI(
    actionButton("filterSCE", "Filter")
  )
)

