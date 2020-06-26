#1GB max upload size
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
options(useFancyQuotes = FALSE)
options(shiny.autoreload = TRUE)

internetConnection <- suppressWarnings(Biobase::testBioCConnection())

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  # library(fs)
  # library(shinyFiles)

  #-----------------------------------------------------------------------------
  # MISC - Used throughout app
  #-----------------------------------------------------------------------------

  #reactive values object
  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset"),
    original = getShinyOption("inputSCEset"),
    batchCorrStatus = "",
    batchResAssay = NULL,
    batchResReddim = NULL,
    diffexgenelist = NULL,
    gsvaRes = NULL,
    gsvaLimma = NULL,
    visplotobject = NULL,
    enrichRes = NULL,
    mastFMRCbin = 0,
    diffexheatmapplot = NULL,
    diffexBmName = NULL,
    celdaMod = NULL,
    celdaList = NULL,
    celdaListAll = NULL,
    celdaListAllNames = NULL,
    celdatSNE = NULL,
    celdaModuleFeature = NULL,
    dimRedPlot = NULL,
    dimRedPlot_geneExp = NULL,
    dendrogram = NULL,
    pcX = NULL,
    pcY = NULL,
    showAssayDetails = FALSE
  )

  #reactive list to store names of results given by the user.
  diffExValues <- reactiveValues(
    index = 0
  )


  #Update all of the columns that depend on pvals columns
  updateColDataNames <- function(){
    pdataOptions <- colnames(colData(vals$counts))
    updateSelectInput(session, "filteredSample",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "deleterowdatacolumn",
                      choices = pdataOptions)
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", pdataOptions))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", pdataOptions))
    updateSelectInput(session, "selectDiffexCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "scMergeCT",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "combatCond",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "batchCorrVar",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "batchCheckVar",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "batchCheckCond",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "mastC1Class",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "mastC2G1Col",
                      choices = c(pdataOptions))
    updateSelectInput(session, "mastC2G2Col",
                      choices = c(pdataOptions))
    updateSelectInput(session, "hurdlecondition",
                      choices = pdataOptions)
    updateSelectInput(session, "mastHMcolData",
                      choices = pdataOptions)
    updateSelectInput(session, "mastHMSplitCol",
                      choices = c('None', 'condition', pdataOptions))
    updateSelectInput(session, "mastFMCluster", choices = pdataOptions)
    updateSelectInput(session, "mastFMHMcolData",
                      choices = pdataOptions)
    updateSelectInput(session, "pathwayPlotVar",
                      choices = pdataOptions)
    updateSelectInput(session, "selectReadDepthCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "selectCellNumCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "selectSnapshotCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "annotModifyChoice",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "visCondn",
                      choices = c("none", pdataOptions))
  }

  updateGeneNames <- function(){
    selectthegenes <- rownames(vals$counts)
    updateSelectizeInput(session, "colorGenes",
                         choices = selectthegenes, server = TRUE)
    updateSelectizeInput(session, "selectvisGenes",
                         choices = selectthegenes, server = TRUE)
    updateSelectizeInput(session, "enrichGenes",
                         choices = selectthegenes, server = TRUE)
  }

  updateFeatureAnnots <- function(){
    selectRowData <- colnames(rowData(vals$counts))
    updateSelectInput(session, "filteredFeature",
                      choices = c("none", selectRowData))
    updateSelectInput(session, "mastHMrowData",
                      choices = selectRowData)
    updateSelectInput(session, "mastHMSplitRow",
                      choices = c('None', 'regulation', selectRowData))
    updateSelectInput(session, "mastFMHMrowData",
                      choices = selectRowData)
  }

  updateNumSamples <- function(){
    numsamples <- ncol(vals$counts)
    updateSelectInput(session, "Knumber",
                      choices = 1:numsamples)
    updateSelectInput(session, "Cnumber",
                      choices = 1:numsamples)
    # updateSelectInput(session, "pcX",
    #                   choices = paste("PC", 1:numsamples, sep = ""),
    #                   selected = "PC1")
    # updateSelectInput(session, "pcY",
    #                   choices = paste("PC", 1:numsamples, sep = ""),
    #                   selected = "PC2")
    updateNumericInput(session, "downsampleNum", value = numsamples,
                       max = numsamples)
  }

  updateAssayInputs <- function(){
    currassays <- names(assays(vals$counts))
    updateSelectInput(session, "dimRedAssaySelect", choices = currassays)
    updateSelectInput(session, "batchCorrAssay", choices = currassays)
    updateSelectInput(session, "batchCheckAssay", choices = currassays)
    updateSelectInput(session, "batchCheckOrigAssay", choices = currassays)
    updateSelectInput(session, "batchCheckCorrAssay",
                      choices = c("", vals$batchResAssay))
    updateSelectInput(session, "diffexAssay", choices = currassays)
    updateSelectInput(session, "mastAssay", choices = currassays)
    updateSelectInput(session, "mastFMAssay", choices = currassays)
    updateSelectInput(session, "mastFMHMAssay", choices = currassays)
    updateSelectInput(session, "pathwayAssay", choices = currassays)
    updateSelectInput(session, "modifyAssaySelect", choices = currassays)
    updateSelectInput(session, "normalizeAssaySelect", choices = currassays)
    updateSelectInput(session, "seuratSelectNormalizationAssay", choices = currassays)
    updateSelectInput(session, "assaySelectFS", choices = currassays)
    updateSelectInput(session, "filterAssaySelect", choices = currassays)
    # updateSelectInput(session, "qcAssaySelect", choices = currassays)
    updateSelectInput(session, "visAssaySelect", choices = currassays)
    updateSelectInput(session, "enrichAssay", choices = currassays)
    updateSelectInput(session, "celdaAssay", choices = currassays)
    updateSelectInput(session, "celdaAssayGS", choices = currassays)
    updateSelectInput(session, "celdaAssaytSNE", choices = currassays)
    updateSelectInput(session, "celdaAssayProbabilityMap",
      choices = currassays)
    updateSelectInput(session, "celdaAssayModuleHeatmap",
      choices = currassays)
    updateSelectInput(session, "depthAssay", choices = currassays)
    updateSelectInput(session, "cellsAssay", choices = currassays)
    updateSelectInput(session, "snapshotAssay", choices = currassays)
  }

  updateReddimInputs <- function(){
    currreddim <- names(reducedDims(vals$counts))
    updateSelectInput(session, "delRedDimType", choices = currreddim)
    updateSelectInput(session, "FastMNNReddim", choices = currreddim)
    updateSelectInput(session, "HarmonyReddim", choices = currreddim)
    updateSelectInput(session, "batchCheckCorrReddim",
                      choices = c("", vals$batchResReddim))
  }

  updateEnrichDB <- function(){
    if (internetConnection){
      enrDB <- enrichR::listEnrichrDbs()$libraryName
    } else {
      enrDB <- ""
    }
    updateSelectInput(session, "enrichDb", choices = c("ALL", enrDB))
  }

  observeEvent(input$consoleToggle, {
    toggle(id = "console")
  })


  # js$disableTabs()

  # Close app on quit
  # session$onSessionEnded(stopApp)

  #-----------------------------------------------------------------------------
  # Page 1: Upload
  #-----------------------------------------------------------------------------

  # Upload data through shiny app

  # Components for uploading directories if user is importing from a preprocessing step
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "base", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFiles::shinyDirChoose(input, "sample", roots = volumes, session = session, restrictions = system.file(package = "base"))

  sample <- reactive(input$sample)
  output$sample <- renderText({
    shinyFiles::parseDirPath(volumes, sample())
  })
  sampleFile <- reactive(input$sampleFile)
  output$sampleFile <- renderText({
    shinyFiles::parseFilePaths(volumes, sampleFile())$datapath
  })
  output$base = renderText({
    shinyDirectoryInput::readDirectoryInput(session, 'directory')
  })
  importCR2Files <- reactiveValues(files = list(), id_count = 0)
  importCR3Files <- reactiveValues(files = list(), id_count = 0)
  importSSFiles <- reactiveValues(files = list(), id_count = 0)
  importBUSFiles <- reactiveValues(files = list(), id_count = 0)
  importSEQFiles <- reactiveValues(files = list(), id_count = 0)
  importOptFiles <- reactiveValues(files = list(), id_count = 0)

  # modal to import all preprocessed data except for CellRanger data
  importModal <- function(failed=FALSE, needsDir=FALSE) {
    modalDialog(
      h3("Sample Name"),
      textInput("sampleName", "*This is the name you would like to give your sample."),
      # only some functions need this input
      if (needsDir)
        h3("Sample ID"),
      if (needsDir)
        textInput("sampleID", "*This name must match your sample's directory name."),

      h3("Base Directory"),
      shinyDirectoryInput::directoryInput('directory', label = 'Choose Directory', value = '~'),
      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("modalOk", "OK")
      )
    )
  }

  # modal to import CellRanger data
  importCRModal <- function() {
    modalDialog(
      h3("Add a Cell Ranger Sample"),
      tags$br(),
      h4("Option 1 - Select a directory containing multiple sample directories (and no other directories)."),
      actionButton("crOpt1", "Add"),
      tags$br(),
      h4("Option 2 - Select a single sample directory."),
      actionButton("crOpt2", "Add"),
      tags$br(),
      h4("Option 3 - Select a directory containing your data files (barcodes.tsv, features.tsv, matrix.mtx)."),
      actionButton("crOpt3", "Add"),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("crOK", "OK")
      )
    )
  }
  # Upload a sample directory (CR) (parent of 'outs' directory)
  importCRSDir <- function(failed = FALSE) {
    modalDialog(
      h3("Sample Directory"),
      shinyDirectoryInput::directoryInput('sDirectory', label = 'Choose Directory', value = '~'),
      h3("Sample Name"),
      h5("If you do not provide an alternate sample name, the sample name will be set to the sample directory name."),
      textInput("sSampleID", ""),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("SDirOK", "OK")
      )
    )
  }
  # Upload a data directory (CR) (parent of 'data files')
  importCRDDir <- function(failed = FALSE) {
    modalDialog(
      h3("Data Directory"),
      shinyDirectoryInput::directoryInput('directory', label = 'Choose Directory', value = '~'),
      h3("Sample Name"),
      textInput("dSampleID", "*This field is mandatory when uploading a data directory"),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("DDirOK", "OK")
      )
    )
  }
  # Upload a base directory (CR) (parent of possibly multiple sample directories)
  importCRBDir <- function(failed = FALSE) {
    modalDialog(
      h3("Base Directory"),
      shinyDirectoryInput::directoryInput('bDirectory', label = 'Choose Directory', value = '~'),
      wellPanel(h5("*For any sample names that you do not provide, the sample name will be set to the sample directory name.")),

      tags$div(id = "bDirTable"),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("BDirOK", "OK")
      )
    )
  }

  # see https://github.com/wleepang/shiny-directory-input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'directory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'directory', value = path)
      }
    }
  )

  # see https://github.com/wleepang/shiny-directory-input
  # for sample directory modal
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$sDirectory
    },
    handlerExpr = {
      if (input$sDirectory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'sDirectory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'sDirectory', value = path)
        if (!is.na(path)) {
          updateTextInput(session, "sSampleID", value = basename(path))
        }
      }
    }
  )

  # event listener for the base directory modal (need to populate table for sample names)
  # see https://github.com/wleepang/shiny-directory-input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$bDirectory
    },
    handlerExpr = {
      if (input$bDirectory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'bDirectory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'bDirectory', value = path)
        # clear the previous table of sample names
        prevPath <- shinyDirectoryInput::readDirectoryInput(session, 'bDirectory')
        count <- 0
        for (prev in list.dirs(prevPath, recursive = FALSE)) {
          count <- count+1
          removeUI(
            selector = paste0("#sampleRow", count)
          )
        }
        # create a new table for the selected directory
        count <- 0
        if (!is.na(path)) {
          counts <- vector()
          for (sample in list.dirs(path, recursive = FALSE)) {
            count <- count+1
            counts <- c(counts, count)
            insertUI(
              selector = "#bDirTable",
              ui = fluidRow(
                id = paste0("sampleRow", count),
                column(6, basename(sample)),
                column(6, textAreaInput(paste0("sampleName", count), "Sample Name", resize = "none", value = basename(sample)))
              )
            )
          }
        }
      }
    }
  )

  # event listeners for "Add Sample" buttons
  observeEvent(input$addCR2Sample, {
    showModal(importCRModal())
  })
  observeEvent(input$crOpt1, {
    removeModal()
    showModal(importCRBDir())
  })
  observeEvent(input$crOpt2, {
    removeModal()
    showModal(importCRSDir())
  })
  observeEvent(input$crOpt3, {
    removeModal()
    showModal(importCRDDir())
  })
  observeEvent(input$addCR3Sample, {
    showModal(importCRModal())
  })
  observeEvent(input$addSSSample, {
    showModal(importModal())
  })
  observeEvent(input$addBUSSample, {
    showModal(importModal())
  })
  observeEvent(input$addSEQSample, {
    showModal(importModal(needsDir = TRUE))
  })
  observeEvent(input$addOptSample, {
    showModal(importModal())
  })

  # function to clear all uploaded files from vectors and UI
  clearAllFiles <- function(fileReactive) {
    for (entry in fileReactive$files) {
      removeUI(selector = paste0("#", entry$id))
    }
    fileReactive$files <- list()
  }

  # event listeners for "Remove Sample" buttons
  observeEvent(input$clearAllCR2, {
    clearAllFiles(importCR2Files)
  })
  observeEvent(input$clearAllCR3, {
    clearAllFiles(importCR3Files)
  })
  observeEvent(input$clearAllSS, {
    clearAllFiles(importSSFiles)
  })
  observeEvent(input$clearAllBUS, {
    clearAllFiles(importBUSFiles)
  })
  observeEvent(input$clearAllSEQ, {
    clearAllFiles(importSEQFiles)
  })
  observeEvent(input$clearAllOpt, {
    clearAllFiles(importOptFiles)
  })

  # event listeners for Cell Ranger import modals' OK buttons
  # sample directory
  observeEvent(input$SDirOK, {
    samplePath <- shinyDirectoryInput::readDirectoryInput(session, 'sDirectory')
    # make sure a directory is selected
    if (identical(samplePath, character(0))) {
      showModal(importCRSDir(failed = TRUE))
    } else {
      # add the files to the appropriate reactiveValues
      if (input$algoChoice == "cellRanger2") {
        id <- paste0("snewSampleCR2", importCR2Files$id_count)
        entry <- list(isDataFile = FALSE, base = paste0(dirname(samplePath), "/"),
                      sample = basename(samplePath), name = input$sSampleID, id = id)
        importCR2Files$files <- c(importCR2Files$files, list(entry))
        importCR2Files$id_count <- importCR2Files$id_count + 1
        selector <- "#newSampleCR2"
      } else {
        id <- paste0("snewSampleCR3", importCR3Files$id_count)
        entry <- list(isDataFile = FALSE, base = paste0(dirname(samplePath), "/"),
                      sample = basename(samplePath), name = input$sSampleID, id = id)
        importCR3Files$files <- c(importCR3Files$files, list(entry))
        importCR3Files$id_count <- importCR3Files$id_count + 1
        selector <- "#newSampleCR3"
      }
      fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
      removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
      # add row for new sample in UI
      insertUI(
        selector = selector,
        ui = fluidRow(
          id = id,
          tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
          column(3, dirname(samplePath)),
          column(3, basename(samplePath)),
          column(3, input$sSampleID),
          column(3, actionButton(paste0("remove", id), "X"))
        )
      )
      # handler to remove the sample that was just added
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        # based on algoChoice, create vector saying which items to keep
        # remove appropriate item from appropriate reactiveValues
        if (input$algoChoice == "cellRanger2") {
          toRemove <- vector()
          for (entry in importCR2Files$files) {
            if (entry$id == id) {
              toRemove <- c(toRemove, FALSE)
            } else {
              toRemove <- c(toRemove, TRUE)
            }
          }
          importCR2Files$files <- importCR2Files$files[toRemove]
        } else {
          toRemove <- vector()
          for (entry in importCR3Files$files) {
            if (entry$id == id) {
              toRemove <- c(toRemove, FALSE)
            } else {
              toRemove <- c(toRemove, TRUE)
            }
          }
          importCR3Files$files <- importCR3Files$files[toRemove]
        }
      })
      removeModal()
    }
  })

  # data directory
  observeEvent(input$DDirOK, {
    dataPath <- shinyDirectoryInput::readDirectoryInput(session, 'directory')
    if ((!nzchar(input$dSampleID)) || (identical(dataPath, character(0)))) {
      showModal(importCRDDir(failed = TRUE))
    } else {
      if (input$algoChoice == "cellRanger2") {
        id <- paste0("dnewSampleCR2", importCR2Files$id_count)
        entry <- list(isDataFile = TRUE, base = "", sample = dataPath, name = input$dSampleID, id = id)
        importCR2Files$files <- c(importCR2Files$files, list(entry))
        importCR2Files$id_count <- importCR2Files$id_count + 1
        selector <- "#newSampleCR2"
      } else {
        id <- paste0("dnewSampleCR3", importCR3Files$id_count)
        entry <- list(isDataFile = TRUE, base = "", sample = dataPath, name = input$dSampleID, id = id)
        importCR3Files$files <- c(importCR3Files$files, list(entry))
        importCR3Files$id_count <- importCR3Files$id_count + 1
        selector <- "#newSampleCR3"
      }
      fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
      removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
      insertUI(
        selector = selector,
        ui = fluidRow(
          id = id,
          tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
          column(3, dataPath),
          column(3, ""),
          column(3, input$dSampleID),
          column(3, actionButton(paste0("remove", id), "X"))
        )
      )
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        if (input$algoChoice == "cellRanger2") {
          toRemove <- vector()
          for (entry in importCR2Files$files) {
            if (entry$id == id) {
              toRemove <- c(toRemove, FALSE)
            } else {
              toRemove <- c(toRemove, TRUE)
            }
          }
          importCR2Files$files <- importCR2Files$files[toRemove]
        } else {
          toRemove <- vector()
          for (entry in importCR3Files$files) {
            if (entry$id == id) {
              toRemove <- c(toRemove, FALSE)
            } else {
              toRemove <- c(toRemove, TRUE)
            }
          }
          importCR3Files$files <- importCR3Files$files[toRemove]
        }
      })
      removeModal()
    }
  })

  # base directory
  observeEvent(input$BDirOK, {
    basePath <- shinyDirectoryInput::readDirectoryInput(session, 'bDirectory')
    # if the user doesn't specify a base directory, show the modal again with the warning message
    if (identical(basePath, character(0))) {
      showModal(importCRBDir(failed = TRUE))
    } else {
      allDirs <- list.dirs(basePath, recursive = FALSE)
      # if we are adding a new CellRangerV2 sample
      if (input$algoChoice == "cellRanger2") {
        allUI <- vector()
        allIDs <- vector()
        count <- 0
        for (sample in allDirs) {
          count <- count + 1
          name <- input[[paste0("sampleName", count)]]
          if (!nzchar(name)) {
            name <- basename(sample)
          }
          id <- paste0("bnewSampleCR2", importCR2Files$id_count)
          entry <- list(isDataFile = FALSE, base = substr(basePath, 1, nchar(basePath)-1), sample = basename(sample), name = name, id = id)
          importCR2Files$files <- c(importCR2Files$files, list(entry))
          fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
          removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
          ui_i <- fluidRow(
            id = id,
            tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
            column(3, basePath),
            column(3, basename(sample)),
            column(3, name),
            column(3, actionButton(paste0("remove", id), "X"))
          )
          importCR2Files$id_count <- importCR2Files$id_count + 1
          allUI <- c(allUI, list(ui_i))
          allIDs <- c(allIDs, id)
        }
        selector <- "#newSampleCR2"
      } else { # if we are adding a new CellRangerV3 sample
        allUI <- vector()
        allIDs <- vector()
        count <- 0
        for (sample in allDirs) {
          count <- count + 1
          name <- input[[paste0("sampleName", count)]]
          if (!nzchar(name)) {
            name <- basename(sample)
          }
          id <- paste0("bnewSampleCR3", importCR3Files$id_count)
          entry <- list(isDataFile = FALSE, base = substr(basePath, 1, nchar(basePath)-1), sample = basename(sample), name = name, id = id)
          importCR3Files$files <- c(importCR3Files$files, list(entry))
          fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
          removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
          ui_i <- fluidRow(
            id = id,
            tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
            column(3, basePath),
            column(3, basename(sample)),
            column(3, name),
            column(3, actionButton(paste0("remove", id), "X"))
          )
          importCR3Files$id_count <- importCR3Files$id_count + 1
          allUI <- c(allUI, list(ui_i))
          allIDs <- c(allIDs, id)
        }
        selector <- "#newSampleCR3"
      }
      # insert all the new sample rows
      for (i in seq_along(allUI)) {
        insertUI(
          selector = selector,
          ui = allUI[i]
        )
      }
      # create event handlers for all the remove buttons
      # from: https://stackoverflow.com/questions/40038749/r-shiny-how-to-write-loop-for-observeevent
      lapply(
        X = allIDs,
        FUN = function(id_i){
          observeEvent(input[[paste0("remove", id_i)]], {
            removeUI(
              selector = paste0("#", id_i)
            )
            if (input$algoChoice == "cellRanger2") {
              toRemove <- vector()
              for (entry in importCR2Files$files) {
                if (entry$id == id_i) {
                  toRemove <- c(toRemove, FALSE)
                } else {
                  toRemove <- c(toRemove, TRUE)
                }
              }
              importCR2Files$files <- importCR2Files$files[toRemove]
            } else {
              toRemove <- vector()
              for (entry in importCR3Files$files) {
                if (entry$id == id_i) {
                  toRemove <- c(toRemove, FALSE)
                } else {
                  toRemove <- c(toRemove, TRUE)
                }
              }
              importCR3Files$files <- importCR3Files$files[toRemove]
            }
          })
        }
      )
      removeModal()
    }
  })

  # event handler for pressing OK on the import modal
  observeEvent(input$modalOk, {
    samplePath <- shinyFiles::parseDirPath(volumes, input$sample)
    basePath <- shinyDirectoryInput::readDirectoryInput(session, 'directory')
    curFiles <- list()
    if ((!nzchar(input$sampleName)) || (identical(basePath, character(0)))) {
      showModal(importModal(failed = TRUE))
    } else {
      entry <- list()
      if (input$algoChoice == "starSolo") {
        curFiles <- importSSFiles
        id <- paste0("newSampleSS", importSSFiles$id_count)
        entry <- list(base = basePath, name = input$sampleName, id = id)
        importSSFiles$files <- c(importSSFiles$files, list(entry))
        importSSFiles$id_count <- importSSFiles$id_count + 1
        selector <- "#newSampleSS"
      } else if (input$algoChoice == "busTools") {
        curFiles <- importBUSFiles
        id <- paste0("newSampleBUS", importBUSFiles$id_count)
        entry <- list(base = basePath, name = input$sampleName, id = id)
        importBUSFiles$files <- c(importBUSFiles$files, list(entry))
        importBUSFiles$id_count <- importBUSFiles$id_count + 1
        selector <- "#newSampleBUS"
      } else if (input$algoChoice == "seqc") {
        curFiles <- importSEQFiles
        id <- paste0("newSampleSEQ", importSEQFiles$id_count)
        entry <- list(base = basePath, sample = input$sampleID, name = input$sampleName, id = id)
        importSEQFiles$files <- c(importSEQFiles$files, list(entry))
        importSEQFiles$id_count <- importSEQFiles$id_count + 1
        selector <- "#newSampleSEQ"
        updateTextInput(session, "sampleID", value = "")
      } else if (input$algoChoice == "optimus") {
        curFiles <- importOptFiles
        id <- paste0("newSampleOpt", importOptFiles$id_count)
        entry <- list(base = basePath, name = input$sampleName, id = id)
        importOptFiles$files <- c(importOptFiles$files, list(entry))
        importOptFiles$id_count <- importOptFiles$id_count + 1
        selector <- "#newSampleOpt"
      }
      fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
      removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
      insertUI(
        selector = selector,
        ui = fluidRow(
          id = id,
          tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
          column(3, basePath),
          column(3, input$sampleID),
          column(3, input$sampleName),
          column(3, actionButton(paste0("remove", id), "X"))
        )
      )
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in curFiles$files) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        curFiles$files <- curFiles$files[toRemove]
      })
      removeModal()
    }
  })

  # Event listener for "Upload" button
  observeEvent(input$uploadData, {
    withBusyIndicatorServer("uploadData", {
      if (input$uploadChoice == "files"){
        vals$original <- importFromFiles(assayFile = input$countsfile$datapath,
                                   annotFile = input$annotFile$datapath,
                                   featureFile = input$featureFile$datapath,
                                   assayName = input$inputAssayType)
      } else if (input$uploadChoice == "example"){
        vals$original <- importExampleData(dataset = input$selectExampleData)
      } else if (input$uploadChoice == "rds") {
        importedrds <- readRDS(input$rdsFile$datapath)
        if (base::inherits(importedrds, "SummarizedExperiment")) {
          vals$original <- importedrds
        } else if (base::inherits(importedrds, "Seurat")) {
          vals$original <- convertSeuratToSCE(importedrds)
        } else {
          showNotification("The '.rds' file should contain a 'SingleCellExperiment' or 'Seurat' object.", type = "error")
        }
      } else if (input$uploadChoice == "directory") {
        # uncomment lines all cbind lines to be able to upload multiple files (and remove the lines above those)
        if (input$algoChoice == "cellRanger2") {
          for (entry in importCR2Files$files) {
            if (entry$isDataFile) {
              sce <- importCellRangerV2Sample(
                sampleDir = entry$sample,
                sampleName = entry$name,
                class = "Matrix",
                delayedArray = FALSE
              )
            } else {
              sce <- importCellRangerV2(
                cellRangerDirs = substr(entry$base, 1, nchar(entry$base)-1),
                sampleDirs = entry$sample,
                sampleNames = entry$name,
                dataType = c("filtered"),
                class = "Matrix",
                delayedArray = FALSE)
            }

            if(is.null(vals$original)) {
              vals$original <- sce
            } else {
              vals$original <- cbind(vals$original, sce)
            }
          }
          clearAllFiles(importCR2Files)
        } else if (input$algoChoice == "cellRanger3") {
          for (entry in importCR3Files$files) {
            if (entry$isDataFile) {
              sce <- importCellRangerV3Sample(
                sampleDir = entry$sample,
                sampleName = entry$name,
                class = "Matrix",
                delayedArray = FALSE
              )
            } else {
              sce <- importCellRangerV3(
                cellRangerDirs = entry$base,
                sampleDirs = entry$sample,
                sampleNames = entry$name,
                dataType = c("filtered"),
                class = "Matrix",
                delayedArray = FALSE)
            }
            if(is.null(vals$original)) {
              vals$original <- sce
            } else {
              vals$original <- cbind(vals$original, sce)
            }
          }
          clearAllFiles(importCR3Files)
        } else if (input$algoChoice == "starSolo") {
          for (entry in importSSFiles$files) {
            sce <- importSTARsolo(
              STARsoloDirs = entry$base,
              samples = entry$name
            )
            vals$original <- sce
            # vals$original <- c(vals$original, list(sce))
          }
          clearAllFiles(importSSFiles)
        } else if (input$algoChoice == "busTools") {
          for (entry in importBUSFiles$files) {
            sce <- importBUStools(
              BUStoolsDirs = substr(entry$base, 1, nchar(entry$base)-1),
              samples = entry$name,
            )
            vals$original <- sce
            # vals$original <- c(vals$original, list(sce))
          }
          clearAllFiles(importBUSFiles)
        } else if (input$algoChoice == "seqc") {
          for (entry in importSEQFiles$files) {
            sce <- importSEQC(
              seqcDirs = entry$base,
              samples = entry$name,
              prefix = entry$sample,
            )
            vals$original <- sce
            # vals$original <- c(vals$original, list(sce))
          }
          clearAllFiles(importSEQFiles)
        } else if (input$algoChoice == "optimus") {
          for (entry in importSEQFiles$files) {
            sce <- importOptimus(
              OptimusDirs = entry$base,
              samples = entry$name
            )
            vals$original <- sce
            # vals$original <- c(vals$original, list(sce))
          }
          clearAllFiles(importOptFiles)
        }
      }

      # Add sample variable if it was not included
      if(is.null(colData(vals$original)$sample)) {
        colData(vals$original)$sample = "sample"
      }

      if (!is.null(vals$original)) {
        vals$counts <- vals$original

        # ToDo: Remove these automatic updates and replace with
        # observeEvents functions that activate upon the tab selection
        updateColDataNames()
        updateNumSamples()
        updateAssayInputs()
        updateGeneNames()
        updateReddimInputs()
        shinyjs::show(id="annotationData")
        js$enableTabs();
      } else {
        shinyalert::shinyalert("Error!", "The data upload failed!",
                               type = "error")
      }
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$diffexheatmapplot <- NULL
      vals$diffexBmName <- NULL
      diffExValues$diffExList <- NULL
      vals$dimRedPlot <- NULL
      vals$dimRedPlot_geneExp <- NULL
      vals$dendrogram <- NULL
      vals$pcX <- NULL
      vals$pcY <- NULL
    })
  })

  #-----------------------------------------------------------------------------
  # Page 2: Data Summary and Filtering
  #-----------------------------------------------------------------------------

  #Sidebar buttons functionality - not an accordion
  shinyjs::onclick("f_hideAllSections", allSections(
    "hide", c(paste("f_collapse", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("f_showAllSections", allSections(
    "show", c(paste("f_collapse", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("f_button1", shinyjs::toggle(id = "f_collapse1",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button2", shinyjs::toggle(id = "f_collapse2",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button3", shinyjs::toggle(id = "f_collapse3",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button4", shinyjs::toggle(id = "f_collapse4",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button5", shinyjs::toggle(id = "f_collapse5",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button6", shinyjs::toggle(id = "f_collapse6",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button7", shinyjs::toggle(id = "f_collapse7",
                                                anim = TRUE), add = TRUE)
  #Button styling
  shinyjs::addClass(id = "f_button1", class = "btn-block")
  shinyjs::addClass(id = "f_button2", class = "btn-block")
  shinyjs::addClass(id = "f_button3", class = "btn-block")
  shinyjs::addClass(id = "f_button4", class = "btn-block")
  shinyjs::addClass(id = "f_button5", class = "btn-block")
  shinyjs::addClass(id = "f_button6", class = "btn-block")
  shinyjs::addClass(id = "f_button7", class = "btn-block")
  shinyjs::addClass(id = "filterData", class = "btn-block")
  shinyjs::addClass(id = "resetData", class = "btn-block")
  shinyjs::addClass(id = "convertGenes", class = "btn-block")
  shinyjs::addClass(id = "deleterowDatabutton", class = "btn-block")
  shinyjs::addClass(id = "downsampleGo", class = "btn-block")

  qc_choice_list <- list("doubletCells", "cxds", "bcds",
                      "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder")

  # Event handler for "Select All" button in QC checklist
  observe({
    if(input$selectallQC == 0) return(NULL)
    else if (input$selectallQC%%2 == 0) {
      updateCheckboxGroupInput(session,"qcAlgos","",choices=qc_choice_list)
    } else {
      updateCheckboxGroupInput(session,"qcAlgos","",choices=qc_choice_list, selected=qc_choice_list)
    }
  })

  qcModal <- function(assays=NULL, geneSetList=FALSE, geneSetListLocation=FALSE,
                      geneSetCollection=FALSE, failed=FALSE, requireAssayStr='') {
    modalDialog(
      h3("QC Paramters - some of the algorithms you have selected require the following extra parameters:"),
      if (!is.null(assays))
        selectInput("qcAssaySelect", paste0("Select assay for ", requireAssayStr), assays),
      if (geneSetList)
        tags$hr(),
      if (geneSetList)
        h4(tags$b("Parameters for QCMetrics:")),
      # The following selectInputs are just place holders until there is gene set code
      if (geneSetList)
        selectInput("geneSetList", "Select Gene Set List", assays),
      if (geneSetListLocation)
        selectInput("geneLocation", "Select Gene Set List Location", assays),
      if (geneSetCollection)
        selectInput("geneCollection", "Select Gene Set Collection", assays),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("modalRunQC", "Run")
      )
    )
  }

  findOverlapping <- function(arr1, arr2) {
    filter <- vector()
    for (x in arr1) {
      if (x %in% arr2) {
        filter <- c(filter, TRUE)
      } else {
        filter <- c(filter, FALSE)
      }
    }
    return(arr1[filter])
  }

  observeEvent(input$runQC, {
    if (is.null(input$qcAlgos)) {
      insertUI(
        selector = "#qcPageErrors",
        ui = wellPanel(id = "noSelected", tags$b("Please select at least one algorithm.", style = "color: red;"))
      )
    } else if (is.null(vals$counts)) {
      insertUI(
        selector = "#qcPageErrors",
        ui = wellPanel(id = "noSCE", tags$b("Please upload a sample first.", style = "color: red;"))
      )
    } else {
      qcAlgosList <- strsplit(input$qcAlgos, " ")
      currassays <- names(assays(vals$counts))
      requireAssay <- list("QCMetrics", "scrublet", "doubletCells", "decontX")
      requireAssayArr <- findOverlapping(qcAlgosList, requireAssay)

      removeUI(
        selector = "#noSelected"
      )
      removeUI(
        selector = "#noSCE"
      )
      if ("QCMetrics" %in% qcAlgosList) {
        showModal(qcModal(assays = currassays, geneSetList = TRUE, geneSetListLocation = TRUE, geneSetCollection = TRUE, requireAssayStr = paste(requireAssayArr, collapse = ', ')))
      } else if (length(requireAssayArr) > 0) {
        showModal(qcModal(assays = currassays, requireAssayStr = paste(requireAssayArr, collapse = ', ')))
      } else {
        runHandler(qcAlgosList)
      }
    }
  })

  observeEvent(input$modalRunQC, {
    qcAlgosList <- strsplit(input$qcAlgos, " ")
    currassays <- names(assays(vals$counts))
    if (is.null(input$qcAssaySelect)) {
      if ("QCMetrics" %in% qcAlgosList) {
        showModal(qcModal(assays = currassays, geneSetList = TRUE, geneSetListLocation = TRUE, geneSetCollection = TRUE, failed= TRUE))
      } else if ("scrublet" %in% qcAlgosList){
        showModal(qcModal(assays = currassays, failed=TRUE))
      } else if ("doubletCells" %in% qcAlgosList) {
        showModal(qcModal(assays = currassays, failed = TRUE))
      } else if ("decontX" %in% qcAlgosList) {
        showModal(qcModal(assays = currassays, failed = TRUE))
      }
    } else {
      removeModal()
      runHandler(qcAlgosList)
    }
  })

  runHandler <- function(qcAlgosList) {
    print(input$qcAssaySelect)
    if ("QCMetrics" %in% qcAlgosList) {
      afterQC <- runCellQC(inSCE = vals$original,
                           algorithms = qcAlgosList,
                           sample = NULL,
                           geneSetList = input$geneSetList,
                           geneSetListLocation = input$geneLocation,
                           geneSetCollection = input$geneCollection,
                           useAssay = input$qcAssaySelect)
    } else if ("scrublet" %in% qcAlgosList){
      afterQC <- runCellQC(inSCE = vals$original,
                           algorithms = qcAlgosList,
                           sample = NULL,
                           useAssay = input$qcAssaySelect)
    } else if ("doubletCells" %in% qcAlgosList) {
      afterQC <- runCellQC(inSCE = vals$original,
                           algorithms = qcAlgosList,
                           sample = NULL,
                           useAssay = input$qcAssaySelect)
    } else if ("decontX" %in% qcAlgosList) {
      afterQC <- runCellQC(inSCE = vals$original,
                           algorithms = qcAlgosList,
                           sample = NULL,
                           useAssay = input$qcAssaySelect)
    } else {
      afterQC <- runCellQC(inSCE = vals$original,
                           algorithms = qcAlgosList,
                           sample = NULL)
    }
    print(afterQC)
    # UNCOMMENT BELOW to show summary table after QC (must uncomment in ui_02_qc as well)
    # output$qcSummary <- renderTable({
    #   req(afterQC)
    #   if(is.null(input$qcAssaySelect)) {
    #     assaySelect <- "counts"
    #   } else {
    #     assaySelect <- input$qcAssaySelect
    #   }
    #   singleCellTK::summarizeTable(inSCE = afterQC,
    #                                useAssay = "counts",
    #                                expressionCutoff = input$minDetectGene)
    # })
    # shinyjs::show(id="qcData")
  }

  #Render data table if there are fewer than 50 samples
  output$contents <- DT::renderDataTable({
    req(vals$counts)
    if (!is.null(getShinyOption("inputSCEset"))){
      updateGeneNames()
    }
    if (!(is.null(vals$counts)) & ncol(vals$counts) < 50){
      temptable <- cbind(rownames(vals$counts), assay(vals$counts, input$filterAssaySelect))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, options = list(scrollX = TRUE), rownames = FALSE)

  #Render histogram of read counts per cell
  output$countshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Arial", size = 14, color = "#7f7f7f")
      x <- list(title = "Reads per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plotly::plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2, function(x) sum(x)),
                      type = "histogram") %>%
        plotly::layout(xaxis = x, yaxis = y)
    } else {
      plotly::plotly_empty(type = "scatter") %>% plotly::add_trace(mode = "lines")
    }
  })

  #Render histogram of genes detected per cell
  output$geneshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Arial", size = 14, color = "#7f7f7f")
      x <- list(title = "Genes detected per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plotly::plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2,
                                function(x) sum(x > 0)), type = "histogram") %>%
        plotly::layout(xaxis = x, yaxis = y)
    } else {
      plotly::plotly_empty(type = "scatter") %>% plotly::add_trace(mode = "lines")
    }
  })

  #random downsample of samples
#  observeEvent(input$downsampleGo, {
#    req(vals$counts)
#    withBusyIndicatorServer("downsampleGo", {
#      vals$counts <- vals$counts[, sample(ncol(vals$counts), input$downsampleNum)]
#      updateNumSamples()
#      vals$diffexheatmapplot <- NULL
#      vals$combatstatus <- ""
#      vals$diffexgenelist <- NULL
#      vals$gsvaRes <- NULL
#      vals$gsvaLimma <- NULL
#      vals$visplotobject <- NULL
#      vals$enrichRes <- NULL
#      vals$diffexBmName <- NULL
#      diffExValues$diffExList <- NULL
#      vals$dimRedPlot <- NULL
#      vals$dimRedPlot_geneExp <- NULL
#      vals$dendrogram <- NULL
#      vals$pcX <- NULL
#      vals$pcY <- NULL
#    })
#  })

  #Render summary table
  output$summarycontents <- renderTable({
    req(vals$counts)

    # Setting 'useAssay=NULL' assumes that the first assay is the one to count
    singleCellTK::summarizeSCE(inSCE = vals$counts,
                               useAssay = NULL,
                               sampleVariableName = "sample")
  }, striped = TRUE, border = TRUE, align = "c", spacing = "l")


  #Reset the data to the original uploaded dataset
  observeEvent(input$resetData, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      vals$counts <- vals$original
      #updateSelectInput(session, "deletesamplelist",
      #                  choices = colnames(vals$counts))
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$diffexBmName <- NULL
      diffExValues$diffExList <- NULL
      vals$dimRedPlot <- NULL
      vals$dimRedPlot_geneExp <- NULL
      vals$dendrogram <- NULL
      vals$pcX <- NULL
      vals$pcY <- NULL
      #Refresh things for the clustering tab
      updateColDataNames()
      updateNumSamples()
      updateAssayInputs()
      updateGeneNames()
      updateEnrichDB()
    }
  })

  #Delete a column from the colData annotations
  observeEvent(input$deleterowDatabutton, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      colData(vals$counts) <- colData(vals$counts)[, !(colnames(colData(vals$counts)) %in% input$deleterowdatacolumn), drop = FALSE]
      updateColDataNames()
    }
  })

  observeEvent(input$filteredSample, {
    output$filterSampleOptions <- renderUI({
      if (input$filteredSample != "none")({
        if (length(unique(colData(vals$counts)[, input$filteredSample])) < 100){
          L <- vector("list", 3)
          L[[1]] <- renderText("Select samples to keep")
          L[[2]] <- wellPanel(style = "overflow-y:scroll; max-height: 100px",
                              list(checkboxGroupInput("filterSampleChoices",
                                                      label = NULL,
                                                      choices = unique(colData(vals$counts)[, input$filteredSample]))),
                              tags$h5(tags$i("Note: the Reset button is in 'Delete Outliers' tab above."))
          )
          L[[3]] <- list(withBusyIndicatorUI(actionButton("runFilterSample", "Filter")))
          return(L)
        } else {
          L <- list(renderText("Annotation must have fewer than 100 options"))
          return(L)
        }
      }) else {
        L <- list()
      }
    })
  })

  #Filter the selected samples
#  observeEvent(input$runFilterSample, {
#    withBusyIndicatorServer("runFilterSample", {
#      filter <- colData(vals$counts)[, input$filteredSample] %in% input$filterSampleChoices
#      vals$counts <- vals$counts[, filter]
#      vals$diffexgenelist <- NULL
#      vals$gsvaRes <- NULL
#      vals$enrichRes <- NULL
#      vals$visplotobject <- NULL
#      vals$diffexheatmapplot <- NULL
#      vals$combatstatus <- ""
#      vals$gsvaLimma <- NULL
#      vals$diffexBmName <- NULL
#      diffExValues$diffExList <- NULL
#      vals$dimRedPlot <- NULL
#      vals$dimRedPlot_geneExp <- NULL
#      vals$dendrogram <- NULL
#      vals$pcX <- NULL
#      vals$pcY <- NULL
#      updateNumSamples()
#    })
#  })

#  observeEvent(input$filteredFeature, {
#    output$filterFeatureOptions <- renderUI({
#      if (input$filteredFeature != "none")({
#        if (length(unique(rowData(vals$counts)[, input$filteredFeature])) < 100){
#          L <- vector("list", 3)
#          L[[1]] <- renderText("Select features to keep")
#          L[[2]] <- wellPanel(style = "overflow-y:scroll; max-height: 100px",
#                              list(checkboxGroupInput("filterFeatureChoices",
#                                                      label = NULL,
#                                                      choices = unique(rowData(vals$counts)[, input$filteredFeature]))))
#          L[[3]] <- list(actionButton("runFilterFeature", "Filter"))
#          return(L)
#        } else {
#          L <- list(renderText("Annotation must have fewer than 100 options"))
#          return(L)
#        }
#      }) else {
#        L <- list()
#      }
#    })
#  })

#  observeEvent(input$orgOrganism, {
#    library(input$orgOrganism, character.only = TRUE)
#    indb <- get(paste(input$orgOrganism))
#    output$orgConvertColumns <- renderUI({
#      tagList(
#        selectInput("orgFromCol", "Select From Annotation:", columns(indb)),
#        selectInput("orgToCol", "Select To Annotation:", columns(indb))
#      )
#    })
#  })

#  observeEvent(input$convertGenes, {
#    req(vals$counts)
#    withBusyIndicatorServer("convertGenes", {
#      vals$counts <- convertGeneIDs(inSCE = vals$counts,
#                                    inSymbol = input$orgFromCol,
#                                    outSymbol = input$orgToCol,
#                                    database = input$orgOrganism)
#      updateGeneNames()
#      vals$diffexgenelist <- NULL
#      vals$gsvaRes <- NULL
#      vals$enrichRes <- NULL
#      vals$visplotobject <- NULL
#      vals$diffexheatmapplot <- NULL
#      vals$diffexBmName <- NULL
#      diffExValues$diffExList <- NULL
#      vals$dimRedPlot <- NULL
#      vals$dimRedPlot_geneExp <- NULL
#      vals$dendrogram <- NULL
#      vals$pcX <- NULL
#      vals$pcY <- NULL
#    })
#  })

  #Filter the selected features
#  observeEvent(input$runFilterFeature, {
#    filter <- rowData(vals$counts)[, input$filteredFeature] %in% input$filterFeatureChoices
#    vals$counts <- vals$counts[filter, ]
#    updateGeneNames()
#    vals$diffexgenelist <- NULL
#    vals$gsvaRes <- NULL
#    vals$enrichRes <- NULL
#    vals$visplotobject <- NULL
#    vals$diffexheatmapplot <- NULL
#    vals$diffexBmName <- NULL
#    diffExValues$diffExList <- NULL
#    vals$dimRedPlot <- NULL
#    vals$dimRedPlot_geneExp <- NULL
#    vals$dendrogram <- NULL
#    vals$pcX <- NULL
#    vals$pcY <- NULL
#  })

  #disable the downloadSCE button if no object is loaded
  isAssayResult <- reactive(is.null(vals$counts))
  observe({
    if (isAssayResult()) {
      shinyjs::disable("downloadSCE")
    } else {
      shinyjs::enable("downloadSCE")
    }
  })

  output$downloadSCE <- downloadHandler(
    filename <- function() {
      paste("SCE-", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
    })

  output$assayList <- renderTable({
    req(vals$counts)
    if (!is.null(vals$counts) & length(names(assays(vals$counts))) > 0){
      data.table(assays = names(assays(vals$counts)))
    }
  })

  output$reducedDimsList <- renderTable({
    req(vals$counts)
    if (!is.null(vals$counts) & length(names(reducedDims(vals$counts))) > 0){
      data.table("Reduced Dimension" = names(reducedDims(vals$counts)))
    }
  })

  observeEvent(input$modifyAssay, {
    req(vals$counts)
    withBusyIndicatorServer("modifyAssay", {
        if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))) {
            showNotification("Assay does not exist!", type = "error")
        } else if (input$modifyAssayOutname == "") {
            showNotification("Assay name cannot be empty!", type = "error")
        } else if (input$modifyAssayOutname %in% names(assays(vals$counts))) {
            showNotification("Assay name already exists! Use another assay name!", type = "error")
        } else {
            if (input$assayModifyAction == "log") {
                if (input$trimAssayCheckbox) {
                    assay(vals$counts, input$modifyAssayOutname) <- log2(assay(vals$counts, input$modifyAssaySelect) + 1)
                    assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
                }
                else {
                    assay(vals$counts, input$modifyAssayOutname) <- log2(assay(vals$counts, input$modifyAssaySelect) + 1)
                }
            }
            else if (input$assayModifyAction == "log1p") {
                if (input$trimAssayCheckbox) {
                    assay(vals$counts, input$modifyAssayOutname) <- log1p(assay(vals$counts, input$modifyAssaySelect))
                    assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
                }
                else {
                    assay(vals$counts, input$modifyAssayOutname) <- log1p(assay(vals$counts, input$modifyAssaySelect))
                }
            }
            else if (input$assayModifyAction == "z.score") {
                if (input$trimAssayCheckbox) {
                    assay(vals$counts, input$modifyAssayOutname) <- computeZScore(assay(vals$counts, input$modifyAssaySelect))
                    assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
                }
                else {
                    assay(vals$counts, input$modifyAssayOutname) <- computeZScore(assay(vals$counts, input$modifyAssaySelect))
                }
            }
            else if(input$assayModifyAction == "trim"){
              assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssaySelect), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
            }
            else {
                showNotification("Error during assay transformation!", type = "error")
            }
          updateAssayInputs()
        }
    })
  })

  observeEvent(input$assayModifyAction,{
    if (input$assayModifyAction == "log"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log"))
    }
    else if (input$assayModifyAction == "log1p"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log1p"))
    }
    else if (input$assayModifyAction == "z.score") {
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Scaled"))
    }
    else if(input$assayModifyAction == "trim"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Trim"))
    }

  })

    observeEvent(input$normalizeAssay, {
    req(vals$counts)
    withBusyIndicatorServer("normalizeAssay", {
      if(!(input$normalizeAssaySelect %in% names(assays(vals$counts)))){
        stop("Selected assay does not exist!")
      }
      else if(input$normalizeAssayOutname == ""){
        stop("Assay Name cannot be empty!")
      }
      else if(input$normalizeAssayOutname %in% names(assays(vals$counts))){
        stop("Your selected Assay Name already exists! Try another Assay Name!")
      }
      else if(input$normalizeAssaySelect == ""){
        stop("Please select an assay before proceeding with normalization!")
      }
      else{
        if (input$normalizeAssayMethodSelect == "LogNormalize"
            || input$normalizeAssayMethodSelect == "CLR"
            || input$normalizeAssayMethodSelect == "RC") {
          vals$counts <- seuratNormalizeData(inSCE = vals$counts,
                                             useAssay = input$normalizeAssaySelect,
                                             normAssayName = input$normalizeAssayOutname,
                                             normalizationMethod = input$normalizeAssayMethodSelect,
                                             scaleFactor = as.numeric(input$normalizationScaleFactor))
          updateAssayInputs()
        }
        else if (input$normalizeAssayMethodSelect == "CPM") {
          assay(vals$counts, input$normalizeAssayOutname) <- scater::calculateCPM(
            x = assay(vals$counts, input$normalizeAssaySelect))
          updateAssayInputs()
        }
        else if(input$normalizeAssayMethodSelect == "LNC"){
          vals$counts <- scater_logNormCounts(
            inSCE = vals$counts,
            logAssayName = input$normalizeAssayOutname,
            useAssay = input$normalizeAssaySelect
          )
          updateAssayInputs()
        }
        else if(input$normalizeAssayMethodSelect == "SCT"){
          vals$counts <- seuratSCTransform(
            inSCE = vals$counts,
            normAssayName = input$normalizeAssayOutname,
            useAssay = input$normalizeAssaySelect
          )
          updateAssayInputs()
        }
      }
    })
})

    observeEvent(input$normalizeAssayMethodSelect, {
      if(input$normalizeAssayMethodSelect == "LogNormalize") {
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratLogNormalize")
      } else if(input$normalizeAssayMethodSelect == "CLR"){
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratCLR")
      } else if(input$normalizeAssayMethodSelect == "RC"){
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratRC")
      } else if(input$normalizeAssayMethodSelect == "CPM"){
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "CPMCounts")
      } else if(input$normalizeAssayMethodSelect == "LNC"){
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "ScaterLogNormCounts")
      } else if(input$normalizeAssayMethodSelect == "SCT"){
        updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratSCTransform")
      }
    })


#  output$colDataDataFrame <- DT::renderDataTable({
#    if (!is.null(vals$counts)){
#      data.frame(colData(vals$counts))
#    }
#  }, options = list(scrollX = TRUE, scrollY = "40vh", pageLength = 30))

  #disable downloadcolData button if the data is not present
  isColDataResult <- reactive(is.null(vals$counts))
  observe({
    if (isColDataResult()) {
      shinyjs::disable("downloadcolData")
    } else {
      shinyjs::enable("downloadcolData")
    }
  })

  #download colData
  output$downloadcolData <- downloadHandler(
    filename = function() {
      paste("colData-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data.frame(colData(vals$counts)), file)
    }
  )

  #upload and replace colData
  observeEvent(input$newAnnotFile, {
    req(input$newAnnotFile)
    req(vals$counts)
    indata <- read.csv(input$newAnnotFile$datapath, header = TRUE, row.names = 1)
    colData(vals$counts) <- DataFrame(indata)
    updateColDataNames()
  })

  #render UI for factor vs numeric
  output$annotModifyUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$annotModifyChoice != "none"){
        HTML(paste("clicking on selected radio button option modifies the selected condition's data in the backend."))
        if (is.factor(colData(vals$counts)[, input$annotModifyChoice])){
          radioButtons("annotTypeSelect", "Field Type:", choices = c("factor", "numeric"), selected = "factor")
        } else {
          radioButtons("annotTypeSelect", "Field Type:", choices = c("factor", "numeric"), selected = "numeric")
        }
      }
    }
  })

  #update factor vs numeric for colData
  observeEvent(input$annotTypeSelect, {
    if (input$annotTypeSelect == "factor" & !is.factor(colData(vals$counts)[, input$annotModifyChoice])){
      colData(vals$counts)[, input$annotModifyChoice] <- as.factor(colData(vals$counts)[, input$annotModifyChoice])
    } else if (input$annotTypeSelect == "numeric" & is.factor(colData(vals$counts)[, input$annotModifyChoice])){
      f <- colData(vals$counts)[, input$annotModifyChoice]
      if (any(is.na(as.numeric(levels(f))[f]))){
        shinyalert::shinyalert("Error!", "Cannot convert to numeric.", type = "error")
      } else {
        colData(vals$counts)[, input$annotModifyChoice] <- as.numeric(levels(f))[f]
      }
    }
  })

  output$annotModifyUIHelpText <- renderUI({
    if (!is.null(vals$counts)){
      if (input$annotModifyChoice != "none"){
        HTML(paste(tags$h5("Note: Clicking on selected radio button option modifies the selected condition's data in the backend.")))
      }
    }
  })

  #-----------------------------------------------------------------------------
  # Page 3: DR & Clustering
  #-----------------------------------------------------------------------------

  #Sidebar buttons functionality - not an accordion
  shinyjs::onclick("c_button2", shinyjs::toggle(id = "c_collapse2",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("c_button3", shinyjs::toggle(id = "c_collapse3",
                                                anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "c_button1", class = "btn-block")
  shinyjs::addClass(id = "c_button2", class = "btn-block")
  shinyjs::addClass(id = "c_button3", class = "btn-block")
  observeEvent(input$delRedDim, {
    req(vals$counts)
    if (!(input$delRedDimType %in% names(reducedDims(vals$counts)))){
      shinyalert::shinyalert("Error!", "reducedDim does not exist!",
                             type = "error")
    } else {
      withBusyIndicatorServer("delRedDim", {
        reducedDim(vals$counts, input$delRedDimType) <- NULL
        updateReddimInputs()
      })
    }
  })

  observeEvent(input$runDimred, {
    if (!is.null(vals$counts)){
      withBusyIndicatorServer("runDimred", {
        if (input$dimRedNameInput == ""){
          shinyalert::shinyalert("Error", "enter a reducedDim name", type = "error")
        } #check for named entered and if its a duplicate
        else if (!is.null(input$dimRedNameInput)){
          if (input$dimRedNameInput %in% names(reducedDims(vals$counts))){
            shinyalert::shinyalert("Error", "Name already exists!", type = "error")
          } else {
            if (input$dimRedPlotMethod == "PCA"){
              if (is.null(reducedDim(vals$counts, input$dimRedNameInput))) {
                vals$counts <- getPCA(inSCE = vals$counts,
                                      useAssay = input$dimRedAssaySelect,
                                      reducedDimName = input$dimRedNameInput)
                updateReddimInputs()
              }
            } else if (input$dimRedPlotMethod == "tSNE"){
              if (is.null(reducedDim(vals$counts, input$dimRedNameInput))) {
                vals$counts <- getTSNE(inSCE = vals$counts,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = input$dimRedNameInput,
                                       perplexity = input$perplexityTSNE,
                                       n_iterations = input$iterTSNE)
                updateReddimInputs()
              }
            } else {
              if (is.null(reducedDim(vals$counts, input$dimRedNameInput))) {
                vals$counts <- getUMAP(inSCE = vals$counts,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = input$dimRedNameInput,
                                       n_neighbors = input$neighborsUMAP,
                                       n_iterations = input$iterUMAP,
                                       alpha = input$alphaUMAP
                )
                updateReddimInputs()
              }
            }
          }
        }
      })
    }
  })

  output$usingReducedDims <- renderUI({
    req(vals$counts)
    selectInput("usingReducedDims", "Select Reduced Dimension Data:", names(reducedDims(vals$counts)))
  })

  output$dimRedAxisSettings <- renderUI({
    req(vals$counts)
    req(input$usingReducedDims)
    if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
      pcComponents <- colnames(reducedDim(vals$counts, input$usingReducedDims))
      pcComponentsSelectedX <- pcComponents[1]
      pcComponentsSelectedY <- pcComponents[2]
      tagList(
        checkboxInput("checkAxis", label = "Modify axes", value = FALSE),
        conditionalPanel(
          condition = "input.checkAxis == true",
          h4("Axis Settings"),
          selectInput("pcX", "X Axis:", pcComponents),
          selectInput("pcY", "Y Axis:", pcComponents, selected = pcComponentsSelectedY)
          )
      )
      }
    })

  observeEvent(input$cUpdatePlot, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }  else {
      withBusyIndicatorServer("cUpdatePlot", {
        if (input$colorBy != "Gene Expression") {
          if (input$axisNames == TRUE) {
            if (input$dimRedAxis1 == "" & input$dimRedAxis2 == "") {
              shinyalert::shinyalert("Error", text = "Enter axis names", type = "error")
            } else {
              comp1 <- input$dimRedAxis1
              comp2 <- input$dimRedAxis2
            }
          } else {
            comp1 <- NULL
            comp2 <- NULL
          }
          #shinyjs doesn't have any visibility functions so have used the following conditions
          if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
            vals$pcX <- input$pcX
            vals$pcY <- input$pcY
          } else {
            vals$pcX <- NULL
            vals$pcY <- NULL
          }
            vals$dimRedPlot <- singleCellTK::plotDimRed(inSCE = vals$counts,
                                                        colorBy = input$colorBy,
                                                        shape = input$shapeBy,
                                                        useAssay = input$dimRedAssaySelect,
                                                        reducedDimName = input$usingReducedDims,
                                                        comp1 = comp1,
                                                        comp2 = comp2,
                                                        pcX = vals$pcX,
                                                        pcY = vals$pcY
            )
          }
      })
    }
    })

  # This code is commented out b/c it was causing a major lag whenever anything else was being
  # updated. Maybe it needs to be changed to observeEvent? - Josh
  # observe({
  #   output$geneExpPlot <- renderPlot({
  #   if (input$colorGeneBy == "Manual Input") {
  #     if (is.null(input$colorGenes)){
  #       ggplot2::ggplot() + ggplot2::theme_bw() +
  #         ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
  #         ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
  #     } else {
  #       if (input$axisNames == TRUE) {
  #         if (input$dimRedAxis1 == "" & input$dimRedAxis2 == "") {
  #           shinyalert::shinyalert("Error", text = "Enter axis names", type = "error")
  #         } else {
  #           comp1 <- input$dimRedAxis1
  #           comp2 <- input$dimRedAxis2
  #         }
  #       } else {
  #         comp1 <- NULL
  #         comp2 <- NULL
  #       }
  #       #shinyjs doesn't have any visibility functions so have used the following conditions
  #       if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
  #         vals$pcX <- input$pcX
  #         vals$pcY <- input$pcY
  #       } else {
  #         vals$pcX <- NULL
  #         vals$pcY <- NULL
  #       }
  #       vals$dimRedPlot_geneExp <- singleCellTK::plotBiomarker(inSCE = vals$counts,
  #                                                              gene = input$colorGenes,
  #                                                              binary = input$colorBinary,
  #                                                              shape = input$shapeBy,
  #                                                              useAssay = input$dimRedAssaySelect,
  #                                                              reducedDimName = input$usingReducedDims,
  #                                                              comp1 = comp1, comp2 = comp2,
  #                                                              x = vals$pcX, y = vals$pcY)
  #       vals$dimRedPlot_geneExp
  #     }
  #   }
  # })
  # })

  output$clusterPlot <- renderPlotly({
    req(vals$dimRedPlot)
    plotly::ggplotly(vals$dimRedPlot)
  })

  observeEvent(input$clusterData, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else if (input$clusterName == "") {
      shinyalert::shinyalert("Error!", "Cluster name required.", type = "error")
    } else {
      withBusyIndicatorServer("clusterData", {
        currdimname <- input$usingReducedDims
        if (input$clusteringAlgorithm == "K-Means"){
          data <- getClusterInputData(vals$counts, input$selectClusterInputData,
                                      useAssay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          koutput <- kmeans(data, input$Knumber)
          name <- input$clusterName
          colData(vals$counts)[, paste(name)] <- factor(koutput$cluster)
          updateColDataNames()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        } else if (input$clusteringAlgorithm == "Clara") {
          data <- getClusterInputData(inSCE = vals$counts,
                                      inputData = input$selectClusterInputData,
                                      useAssay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          coutput <- cluster::clara(data, input$Cnumber)
          name <- input$clusterName
          colData(vals$counts)[, paste(name)] <- factor(coutput$clustering)
          updateColDataNames()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        }
      })
    }
  })

  #TODO: this doesn't work with multiple pca dims
  output$pctable <- renderTable({
      if (!is.null(vals$counts)){
       # HTML(tags$h4("PC Table:"))
          if (any(grepl(pattern = "PC*", x = colnames(reducedDim(vals$counts, input$usingReducedDims))))) {
            if (nrow(pcaVariances(vals$counts)) == ncol(vals$counts)){
              data.frame(PC = paste("PC", seq_len(ncol(vals$counts)), sep = ""),
                         Variances = pcaVariances(vals$counts)$percentVar * 100)[1:10, ]
            }
          }
        }
    })

  #Gene visualization
  output$visOptions <- renderUI({
    if (!is.null(vals$counts)){
      if (input$visPlotMethod != "heatmap") {
        tagList(
          checkboxInput("visFWrap", "Plot genes individually?", value = TRUE)
        )
      } else {
        tagList(
          checkboxInput("visScaleHMap", "Scale expression values?", value = TRUE)
        )
      }
    }
  })

  output$visBioGenes <- renderUI({
    if (!is.null(vals$counts)) {
      selectInput("selVisBioGenes", "Select Gene List(s):",
                  names(which(apply(rowData(vals$counts), 2, function(a) length(unique(a)) == 2) == TRUE)))
    }
  })

  observeEvent(input$plotvis, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("plotvis", {
        tryCatch({
          if (input$visCondn == "none"){
            incondition <- NULL
          } else {
            incondition <- input$visCondn
          }
          if (input$visGeneList == "selVisRadioGenes"){
            visGList <- input$selectvisGenes
          } else  {
            visGList <- rownames(vals$counts)[SingleCellExperiment::rowData(vals$counts)[, input$selVisBioGenes] == 1]
            #if Gene list > 25 choose the top 25 which is ordered according to p-val
            if (length(visGList) > 25) {
              visGList <- visGList[1:25]
            }
          }
          vals$visplotobject <- visPlot(inSCE = vals$counts,
                                        useAssay = input$visAssaySelect,
                                        method =  input$visPlotMethod,
                                        condition = incondition,
                                        glist =  visGList,
                                        facetWrap = input$visFWrap,
                                        scaleHMap = input$visScaleHMap)
        },
        error = function(e){
          shinyalert::shinyalert("Error!", e$message, type = "error")
        })
      })
    }
  })

  output$visPlot <- renderPlot({
    req(vals$visplotobject)
    vals$visplotobject
  }, height = 600)

  #-----------------------------------------------------------------------------
  # Page 3.2: Celda
  #-----------------------------------------------------------------------------

  # onclick buttons
  shinyjs::onclick("celdaBasicSet",
                   shinyjs::toggle(id = "celdaCollapse1",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("celdaAdvSet",
                   shinyjs::toggle(id = "celdaCollapse2",
                                   anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "celdaBasicSet", class = "btn-block")
  shinyjs::addClass(id = "celdaAdvSet", class = "btn-block")

  shinyjs::onclick("celdaBasicSetGS",
    shinyjs::toggle(id = "celdaCollapseGS1",
      anim = TRUE), add = TRUE)
  shinyjs::onclick("celdaAdvSetGS",
    shinyjs::toggle(id = "celdaCollapseGS2",
      anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "celdaBasicSetGS", class = "btn-block")
  shinyjs::addClass(id = "celdaAdvSetGS", class = "btn-block")

  # shinyjs::onclick("celdatSNESet",
  #   shinyjs::toggle(id = "celdaCollapsetSNE",
  #     anim = TRUE), add = TRUE)
  # shinyjs::addClass(id = "celdatSNESet", class = "btn-block")


  # celda clustering tab
  observeEvent(input$runCelda, {
    # is there an error or not
    if (is.null(vals$counts)) {
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      # selected count matrix
      cm <- assay(vals$counts, input$celdaAssay)
    }

    # And each row/column of the count matrix must have at least one count
    if (sum(rowSums(cm) == 0) >= 1 | sum(colSums(cm) == 0) >= 1) {
      shinyalert::shinyalert("Error!",
        "Each row and column of the count matrix must have at least one count.
        Filter the data first.",
        type = "error")
    }

    # Ensure that number of genes / cells is never more than
    # the number of requested clusters for each
    if (!is.null(input$cellClusterC) && ncol(cm) < input$cellClusterC) {
      shinyalert::shinyalert("Error!",
        "Number of cells (columns) in count matrix must be >= K",
        type = "error")
    }

    if (!is.null(input$geneModuleG) && nrow(cm) < input$geneModuleG) {
      shinyalert::shinyalert("Error!",
        "Number of genes (rows) in count matrix must be >= L",
        type = "error")
    }

    if (!is.null(input$geneModuleCG) && nrow(cm) < input$geneModuleCG) {
      shinyalert::shinyalert("Error!",
        "Number of genes (rows) in count matrix must be >= L",
        type = "error")
    }

    if (!is.null(input$cellClusterCG) && ncol(cm) < input$cellClusterCG) {
      shinyalert::shinyalert("Error!",
        "Number of cells (columns) in count matrix must be >= K",
        type = "error")
    } else {

      withBusyIndicatorServer("runCelda", {
        if (input$celdaModel == "celda_C") {
          vals$celdaMod <- celda_C(counts = assay(vals$counts,
            input$celdaAssay),
            K = input$cellClusterC,
            alpha = input$celdaAlpha,
            beta = input$celdaBeta,
            algorithm = input$celdaAlgorithm,
            stopIter = input$celdaStopIter,
            maxIter = input$celdaMaxIter,
            splitOnIter = input$celdaSplitIter,
            seed = input$celdaSeed,
            nchains = input$celdaNChains)
          #cores = input$celdaCores)
          colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
          updateColDataNames()

        } else if (input$celdaModel == "celda_G") {
          vals$celdaMod <- celda_G(counts = assay(vals$counts,
            input$celdaAssay),
            L = input$geneModuleG,
            beta = input$celdaBeta,
            delta = input$celdaDelta,
            gamma = input$celdaGamma,
            stopIter = input$celdaStopIter,
            maxIter = input$celdaMaxIter,
            splitOnIter = input$celdaSplitIter,
            seed = input$celdaSeed,
            nchains = input$celdaNChains)
          #cores = input$celdaCores)
          rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
          updateFeatureAnnots()
          # update feature modules in module heatmap panel
          updateSelectInput(session, "celdaFeatureModule",
            choices = seq_len(vals$celdaMod@params$L))

        } else if (input$celdaModel == "celda_CG") {
          vals$celdaMod <- celda_CG(counts = assay(vals$counts,input$celdaAssay),
            K = input$cellClusterCG,
            L = input$geneModuleCG,
            alpha = input$celdaAlpha,
            beta = input$celdaBeta,
            delta = input$celdaDelta,
            gamma = input$celdaGamma,
            algorithm = input$celdaAlgorithm,
            stopIter = input$celdaStopIter,
            maxIter = input$celdaMaxIter,
            splitOnIter = input$celdaSplitIter,
            seed = input$celdaSeed,
            nchains = input$celdaNChains)
          #cores = input$celdaCores)
          colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
          rowData(vals$counts)$celdaGeneModule <- rep(NA, nrow(vals$counts))
          rowData(vals$counts)$celdaGeneModule[ix] <- vals$celdaMod@clusters$y
          updateColDataNames()
          updateFeatureAnnots()
          # update feature modules in module heatmap panel
          updateSelectInput(session, "celdaFeatureModule",
            choices = seq_len(vals$celdaMod@params$L))
          updateSelectInput(session, "celdatSNEFeature",
            choices = vals$celdaMod@names$row)
        }
      })
    }
  })

  # disable the renderCeldaHeatmap button if no celda model is present
  isCeldaModelPresent <- reactive(is.null(vals$celdaMod))
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderHeatmap")
    } else {
      shinyjs::enable("renderHeatmap")
    }
  })


  # show celda heatmap
  observeEvent(input$renderHeatmap, {
    #celdaHeatmap <- eventReactive(input$showHeatmap, {

    # is there an error or not
    if (is.null(vals$celdaMod)) {
      shinyalert::shinyalert("Error!",
        "Celda Model Not Found.",
        type = "error")
    } else {

      withBusyIndicatorServer("renderHeatmap",
        output$celdaHeatmap <- renderPlot({
          g <- celdaHeatmap(counts = assay(vals$counts,input$celdaAssay),
            celdaMod = vals$celdaMod)
          g
        }, height = 600)
      )}
  })


  #disable the downloadSCECelda button if no object is loaded
  isAssayResultCelda <- reactive(is.null(vals$counts))
  observe({
    if (isAssayResultCelda()) {
      shinyjs::disable("downloadSCECelda")
    } else {
      shinyjs::enable("downloadSCECelda")
    }
  })

  output$downloadSCECelda <- downloadHandler(
    filename <- function() {
      paste("SCE_", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
    })

  observeEvent(input$celdaSelectGSList, {
    updateSelectInput(session, "celdaSelectGSMod",
      choices = vals$celdaListAllNames[[input$celdaSelectGSList]])
  })


  # show celda perplexity plot
  observeEvent(input$renderPerplexityPlot, {
    # is there an error or not
    if (is.null(vals$celdaList)) {
      shinyalert::shinyalert("Error!",
        "Celda List Not Found.",
        type = "error")
    } else {
      withBusyIndicatorServer("renderPerplexityPlot", {
        vals$celdaList <- resamplePerplexity(counts = assay(vals$counts,
          input$celdaAssayGS),
          celdaList = vals$celdaList)
        output$celdaPerplexityPlot <- renderPlot({
          g <- plotGridSearchPerplexity(celdaList = vals$celdaList)
          g
        }, height = 600)
      })}
  })

  # Confirm celda model
  # disable the Confirm Selection button if no celda list is present
  isCeldaModelListPresent <- reactive(is.null(vals$celdaListAll))
  observe({
    if (isCeldaModelListPresent()) {
      shinyjs::disable("confirmCeldaModel")
    } else {
      shinyjs::enable("confirmCeldaModel")
    }
  })


  observeEvent(input$confirmCeldaModel, {
    withBusyIndicatorServer("confirmCeldaModel", {
      vals$celdaMod <- vals$celdaListAll[[
        input$celdaSelectGSList]][[input$celdaSelectGSMod]]
      # update data annotations
      if (paste(strsplit(input$celdaSelectGSMod,
        "_")[[1]][seq(2)], collapse  = "_") == "celda_C") {
        colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
        updateColDataNames()
      } else if (paste(strsplit(input$celdaSelectGSMod,
        "_")[[1]][seq(2)], collapse  = "_") == "celda_G") {
        rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
        updateFeatureAnnots()
      } else if (paste(strsplit(input$celdaSelectGSMod,
        "_")[[1]][seq(2)], collapse  = "_") == "celda_CG") {
        colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
        rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
        updateColDataNames()
        updateFeatureAnnots()
        updateSelectInput(session, "celdatSNEFeature",
          choices = vals$celdaMod@names$row)
      }
      # update feature modules in module heatmap panel
      updateSelectInput(session, "celdaFeatureModule",
        choices = seq_len(vals$celdaMod@params$L))
    })
  })

  # Visualize tab
  # tSNE sub-panel
  # disable the runCeldatSNE button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("runCeldatSNE")
    } else {
      shinyjs::enable("runCeldatSNE")
    }
  })

  # run celda tSNE
  observeEvent(input$runCeldatSNE, {
    withBusyIndicatorServer("runCeldatSNE", {
      vals$celdatSNE <- celdaTsne(counts = assay(vals$counts,
        input$celdaAssaytSNE),
        celdaMod = vals$celdaMod,
        maxCells = input$celdatSNEmaxCells,
        minClusterSize = input$celdatSNEminClusterSize,
        perplexity = input$celdatSNEPerplexity,
        maxIter = input$celdatSNEmaxIter,
        seed = input$celdatSNESeed)
    })
  })

  # disable the renderCeldatSNE button if no celda tSNE result is present
  isCeldatSNEResult <- reactive(is.null(vals$celdatSNE))
  observe({
    if (isCeldatSNEResult()) {
      shinyjs::disable("renderCeldatSNEByCellCluster")
      shinyjs::disable("renderCeldatSNEModule")
      shinyjs::disable("renderCeldatSNEFeature")
    } else {
      shinyjs::enable("renderCeldatSNEByCellCluster")
      shinyjs::enable("renderCeldatSNEModule")
      shinyjs::enable("renderCeldatSNEFeature")
    }
  })

  # show celda tSNE colored by cell cluster
  observeEvent(input$renderCeldatSNEByCellCluster, {
      withBusyIndicatorServer("renderCeldatSNEByCellCluster",
        output$celdatSNECellClusterPlot <- renderPlot({
          g <- plotDimReduceCluster(
            dim1 = vals$celdatSNE[, 1],
            dim2 = vals$celdatSNE[, 2],
            cluster = clusters(vals$celdaMod)$z)
          g}, height = 600)
      )}
  )

  # show celda tSNE colored by module probabilities
  observeEvent(input$renderCeldatSNEModule, {
    withBusyIndicatorServer("renderCeldatSNEModule",
      output$celdatSNEModulePlot <- renderPlot({
        g <- plotDimReduceModule(
          dim1 = vals$celdatSNE[, 1],
          dim2 = vals$celdatSNE[, 2],
          counts = assay(vals$counts, input$celdaAssaytSNE),
          celdaMod = vals$celdaMod)
        g}, height = 600)
    )}
  )

  # show celda tSNE colored by feature expression
  observeEvent(input$renderCeldatSNEFeature, {
    withBusyIndicatorServer("renderCeldatSNEFeature",
      output$celdatSNEFeaturePlot <- renderPlot({
        g <- plotDimReduceFeature(
          dim1 = vals$celdatSNE[, 1],
          dim2 = vals$celdatSNE[, 2],
          counts = assay(vals$counts, input$celdaAssaytSNE),
          features = input$celdatSNEFeature)
        g}, height = 600)
    )}
  )

  # Probability Map sub-panel
  # disable the renderCeldaProbabilityMap button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderCeldaProbabilityMap")
    } else {
      shinyjs::enable("renderCeldaProbabilityMap")
    }
  })

  # show celda Probability Map
  observeEvent(input$renderCeldaProbabilityMap, {
    withBusyIndicatorServer("renderCeldaProbabilityMap",
      output$celdaProbabilityMapPlot <- renderPlot({
        g <- celdaProbabilityMap(counts = assay(vals$counts,
          input$celdaAssayProbabilityMap),
          celdaMod = vals$celdaMod)
        g}, height = 600)
    )}
  )


  # Module Heatmap sub-panel
  # disable the renderCeldaModuleHeatmap button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderCeldaModuleHeatmap")
    } else {
      shinyjs::enable("renderCeldaModuleHeatmap")
    }
  })

  # show celda Module Heatmap
  observeEvent(input$renderCeldaModuleHeatmap, {
    withBusyIndicatorServer("renderCeldaModuleHeatmap",
      output$celdaModuleHeatmapPlot <- renderPlot({
        g <- moduleHeatmap(counts = assay(vals$counts,
          input$celdaAssayModuleHeatmap),
          celdaMod = vals$celdaMod,
          featureModule = as.integer(input$celdaFeatureModule),
          topCells = input$celdaModuleTopCells,
          #normalize = as.logical(input$celdaFeatureModuleNormalize),
          showFeaturenames = as.logical(input$celdaModuleFeatureNames))
        g}, height = 600)
    )}
  )


  # download all celda lists
  # disable the downloadCeldaList button if no celda list is present
  isAssayResultCeldaList <- reactive(is.null(vals$celdaListAll))
  observe({
    if (isAssayResultCeldaList()) {
      shinyjs::disable("downloadAllCeldaLists")
    } else {
      shinyjs::enable("downloadAllCeldaLists")
    }
  })

  output$downloadAllCeldaLists <- downloadHandler(
    filename <- function() {
      paste("All_Celda_Lists_", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$celdaListAll, file)
    })


  #-----------------------------------------------------------------------------
  # Page 3.3: Cell Viewer
  #-----------------------------------------------------------------------------
  #-+-+-+-+-+-For Functional Panel collapse##############
  shinyjs::onclick("cv_button1", shinyjs::toggle(id = "cv_collapse1",
    anim = TRUE), add = TRUE)
  shinyjs::onclick("cv_button2", shinyjs::toggle(id = "cv_collapse2",
    anim = TRUE), add = TRUE)
  shinyjs::onclick("cv_button3", shinyjs::toggle(id = "cv_collapse3",
    anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "cv_button1", class = "btn-block")
  shinyjs::addClass(id = "cv_button2", class = "btn-block")
  shinyjs::addClass(id = "cv_button3", class = "btn-block")
  colorbrewer_list <- rownames(RColorBrewer::brewer.pal.info)
  color_table <- RColorBrewer::brewer.pal.info %>% data.frame()
  color_seqdiv <- rownames(color_table[which(color_table$category == "div"
    |color_table$category == "seq"),])

  # #-+-+-+-+-+-For Input Observe##############
  # observe({
  #   # is there an error or not
  #   if (is.null(vals$counts)) {
  #     # shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
  #   } else {
  #     #colorbrewer_list <- rownames(RColorBrewer::brewer.pal.info)
  #     #color_table <- RColorBrewer::brewer.pal.info %>% data.frame()
  #     #color_seqdiv <- rownames(color_table[which(color_table$category == "div"
  #     #                                           |color_table$category == "seq"),])
  #     #from sce
  #     cell_list <- BiocGenerics::colnames(vals$counts)
  #     gene_list <- BiocGenerics::rownames(vals$counts)
  #     #from assays
  #     method_list <- names(assays(vals$counts))
  #     #from reduced
  #     approach_list <- names(reducedDims(vals$counts))
  #     #from colData
  #     annotation_list <- names(colData(vals$counts))
  #
  #     updateSelectInput(session, "QuickAccess",
  #       choices = c("",approach_list,"Custom"))
  #     updateSelectInput(session, "ApproachSelect_Xaxis",
  #       choices = c(approach_list))
  #     updateSelectInput(session, "AdvancedMethodSelect_Xaxis",
  #       choices = c(method_list))
  #     updateSelectInput(session, "GeneSelect_Assays_Xaxis",
  #       choices = c(gene_list))
  #     updateSelectInput(session, "AnnotationSelect_Xaxis",
  #       choices = c(annotation_list))
  #     updateSelectInput(session, "ApproachSelect_Yaxis",
  #       choices = c(approach_list))
  #     updateSelectInput(session, "AdvancedMethodSelect_Yaxis",
  #       choices = c(method_list))
  #     updateSelectInput(session, "GeneSelect_Assays_Yaxis",
  #       choices = c(gene_list))
  #     updateSelectInput(session, "AnnotationSelect_Yaxis",
  #       choices = c(annotation_list))
  #     updateSelectInput(session, "ApproachSelect_Colorby",
  #       choices = c(approach_list))
  #     updateSelectInput(session, "AdvancedMethodSelect_Colorby",
  #       choices = c(method_list))
  #     updateSelectInput(session, "GeneSelect_Assays_Colorby",
  #       choices = c(gene_list))
  #     updateSelectInput(session, "AnnotationSelect_Colorby",
  #       choices = c(annotation_list))
  #     updateSelectizeInput(session, "adjustgroupby", label = NULL, choices = c("None", annotation_list))
  #     updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:",
  #       choices = c("RdYlBu",color_seqdiv))
  #   }
  # })
  #
  # #-+-+-+-+-+-For Advanced Input Observe##############
  # ###ApproachSelect to DimensionSelect X-Axis
  # observe({
  #   if (!is.null(vals$counts)){
  #     len <- length(SingleCellExperiment::reducedDims(vals$counts))
  #     if (!is.null(input$ApproachSelect_Xaxis) & len > 0){
  #       Df <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Xaxis))
  #       xs <- colnames(Df)
  #       updateSelectInput(session, "ColumnSelect_Xaxis", choices = c(xs))
  #       rm(Df)
  #     }
  #   }
  # })
  # ###ApproachSelect to DimensionSelect Y-Axis
  # observe({
  #   if (!is.null(vals$counts)){
  #     len <- length(SingleCellExperiment::reducedDims(vals$counts))
  #     if (!is.null(input$ApproachSelect_Yaxis) & len > 0){
  #       Df2 <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Yaxis))
  #       xs2 <- colnames(Df2)
  #       xs2 <- sort(xs2, decreasing = TRUE)
  #       updateSelectInput(session, "ColumnSelect_Yaxis", choices = c(xs2))
  #       rm(Df2)
  #     }
  #   }
  # })
  # ###ApproachSelect to DimensionSelect Colorby
  # observe({
  #   if (!is.null(vals$counts)){
  #     len <- length(SingleCellExperiment::reducedDims(vals$counts))
  #     if (!is.null(input$ApproachSelect_Colorby) & len > 0){
  #       Df3 <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Colorby))
  #       xs3 <- colnames(Df3)
  #       updateSelectInput(session, "ColumnSelect_Colorby", choices = c(xs3))
  #       rm(Df3)
  #     }
  #   }
  # })
  #
  # #-+-+-+-+-+-Observe Group by###################################################
  # ###Observe Radio Button Select Value Type
  # observe({
  #   if (!is.null(vals$counts)){
  #     if (input$adjustgroupby !=  'None'){
  #       #Integer,level>25#
  #       if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
  #         & length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))>25){
  #         updateRadioButtons(session, "SelectValueType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Continuous")
  #         shinyjs::delay(5,shinyjs::disable("SelectValueType"))
  #         #Integer,level<25#
  #       }else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
  #         & length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))<=25){
  #         updateRadioButtons(session, "SelectValueType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Categorical")
  #         shinyjs::enable("SelectValueType")
  #         #Numeric,noninteger#
  #       }else if(is.numeric(colData(vals$counts)@listData[[input$adjustgroupby]])){
  #         updateRadioButtons(session, "SelectValueType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Continuous")
  #         shinyjs::delay(5,shinyjs::disable("SelectValueType"))
  #         #Categorical#
  #       }else{
  #         updateRadioButtons(session, "SelectValueType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Categorical")
  #         shinyjs::delay(5,shinyjs::disable("SelectValueType"))}
  #     }
  #   }
  # })#observe_end
  #
  # ###Observe Check Box Check Binning & Text Input Number of Bins:
  #
  # observe({
  #   if (!is.null(vals$counts)){
  #     if (input$adjustgroupby !=  'None'){
  #       #Integer,level>25#
  #       if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
  #         &length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))>25){
  #         updateCheckboxInput(session,"checkbinning","Perform Binning", value = TRUE)
  #         shinyjs::delay(5,shinyjs::disable("checkbinning"))
  #         shinyjs::enable("adjustbinning")
  #         #Integer,level<25,continuous
  #       }else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
  #         &length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))<=25
  #         &input$SelectValueType == "Continuous"){
  #         updateCheckboxInput(session,"checkbinning","Perform Binning", value = TRUE)
  #         shinyjs::delay(5,shinyjs::disable("checkbinning"))
  #         shinyjs::enable("adjustbinning")
  #         #Integer,level<25,Categorical
  #       }else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
  #         &length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))<=25
  #         &input$SelectValueType == "Categorical"){
  #         updateCheckboxInput(session,"checkbinning","Perform Binning", value = FALSE)
  #         shinyjs::delay(5,shinyjs::disable("checkbinning"))
  #         shinyjs::disable("adjustbinning")
  #         #Numeric,noninteger
  #       }else if(is.numeric(colData(vals$counts)@listData[[input$adjustgroupby]])){
  #         updateCheckboxInput(session,"checkbinning","Perform Binning", value = TRUE)
  #         shinyjs::delay(5,shinyjs::disable("checkbinning"))
  #         shinyjs::enable("adjustbinning")
  #         #Categorical
  #       }else{updateCheckboxInput(session,"checkbinning","Perform Binning", value = FALSE)
  #         shinyjs::delay(5,shinyjs::disable("checkbinning"))
  #         shinyjs::disable("adjustbinning")
  #       }
  #     }
  #   }
  # })#observe_end
  #
  # #-+-+-+-+-+-Observe Color bye###################################################
  # ###Observe Radio Button Select Value Type
  # observe({
  #   if (!is.null(vals$counts)){
  #     if (input$TypeSelect_Colorby != 'Pick a Color'){
  #       ###If Cell Annotation###############################################################
  #       if(input$TypeSelect_Colorby == 'Cell Annotation'){
  #         ###If Cell Annotation numeric
  #         if(!is.numeric(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])){
  #           updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #             choices = c("Categorical", "Continuous"),
  #             selected = "Categorical")
  #           shinyjs::delay(5,shinyjs::disable("SelectColorType"))
  #
  #
  #         }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
  #           &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25){
  #           updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #             choices = c("Categorical", "Continuous"),
  #             selected = "Categorical")
  #           shinyjs::enable("SelectColorType")
  #
  #         }else{updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Continuous")
  #           shinyjs::delay(5,shinyjs::disable("SelectColorType"))}
  #
  #         ###If ReducedData##########################################################
  #       }else if(input$TypeSelect_Colorby == 'Reduced Dimensions'){
  #         Dfcolor <- data.frame(reducedDims(vals$counts)@listData[[input$ApproachSelect_Colorby]])
  #         if(input$ColumnSelect_Colorby %in% colnames(Dfcolor)){
  #           Dfcolor <- Dfcolor[which(colnames(Dfcolor) == input$ColumnSelect_Colorby)]
  #           ###If ReducedData numeric
  #
  #           if(!is.numeric(Dfcolor[,1])){
  #             updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #               choices = c("Categorical", "Continuous"),
  #               selected = "Categorical")
  #             shinyjs::delay(5,shinyjs::disable("SelectColorType"))
  #
  #
  #           }else if(is.integer(Dfcolor[,1])
  #             &length(levels(as.factor(Dfcolor[,1])))<=25){
  #             updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #               choices = c("Categorical", "Continuous"),
  #               selected = "Categorical")
  #             shinyjs::enable("SelectColorType")
  #
  #           }else{updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #             choices = c("Categorical", "Continuous"),
  #             selected = "Continuous")
  #             shinyjs::delay(5,shinyjs::disable("SelectColorType"))}
  #         }
  #         ###If Expression Assays###########################################################
  #       }else{Dfassay <- assay(vals$counts, input$AdvancedMethodSelect_Colorby)
  #       if(input$GeneSelect_Assays_Colorby %in% rownames(Dfassay)){
  #         Dfassay <- data.frame(Dfassay[which(rownames(Dfassay)== input$GeneSelect_Assays_Colorby),])
  #
  #         if(!is.numeric(Dfassay[,1])){
  #           updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #             choices = c("Categorical", "Continuous"),
  #             selected = "Categorical")
  #           shinyjs::delay(5,shinyjs::disable("SelectColorType"))
  #
  #
  #         }else if(is.integer(Dfassay[,1])
  #           &length(levels(as.factor(Dfassay[,1])))<=25){
  #           updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #             choices = c("Categorical", "Continuous"),
  #             selected = "Categorical")
  #           shinyjs::enable("SelectColorType")
  #
  #         }else{updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
  #           choices = c("Categorical", "Continuous"),
  #           selected = "Continuous")
  #           shinyjs::delay(5,shinyjs::disable("SelectColorType"))}
  #       }
  #       }
  #     }
  #   }
  # })###observe_end
  #
  # ###Observe Check Box Check Binning & Text Input Number of Bins:
  # observe({
  #   if (!is.null(vals$counts)){
  #     ###If Cell Annotation###############################################################
  #     if(input$TypeSelect_Colorby != 'Pick a Color'){
  #
  #       if(input$TypeSelect_Colorby == 'Cell Annotation'){
  #         if(!is.numeric(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])){
  #           updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #           shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #           shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #           updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #         }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
  #           &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25
  #           &input$SelectColorType == 'Categorical'){
  #           updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #           shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #           shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #           updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #         }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
  #           &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25
  #           &input$SelectColorType == 'Continuous'){
  #
  #           shinyjs::enable("checkColorbinning")
  #           if(input$checkColorbinning == TRUE){
  #             shinyjs::enable("adjustColorbinning")
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #           else{
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))}
  #
  #         }else{
  #
  #           shinyjs::enable("checkColorbinning")
  #           if(input$checkColorbinning == TRUE){
  #             shinyjs::enable("adjustColorbinning")
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #           else{
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))}
  #         }
  #
  #
  #         ###If Reduce Dimensions##############################################################
  #       }else if(input$TypeSelect_Colorby == 'Reduced Dimensions'){
  #         Dfcolor <- data.frame(reducedDims(vals$counts)@listData[[input$ApproachSelect_Colorby]])
  #         if(input$ColumnSelect_Colorby %in% colnames(Dfcolor)){
  #           Dfcolor <- Dfcolor[which(colnames(Dfcolor) == input$ColumnSelect_Colorby)]
  #
  #           if(!is.numeric(Dfcolor[,1])){
  #             updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #             shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #           }else if(is.integer(Dfcolor[,1])
  #             &length(levels(as.factor(Dfcolor[,1])))<=25
  #             &input$SelectColorType == 'Categorical'){
  #             updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #             shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #           }else if(is.integer(Dfcolor[,1])
  #             &length(levels(as.factor(Dfcolor[,1])))<=25
  #             &input$SelectColorType == 'Continuous'){
  #
  #             shinyjs::enable("checkColorbinning")
  #             if(input$checkColorbinning == TRUE){
  #               shinyjs::enable("adjustColorbinning")
  #               updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #             else{
  #               shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #               updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))}
  #
  #           }else{
  #
  #             shinyjs::enable("checkColorbinning")
  #             if(input$checkColorbinning == TRUE){
  #               shinyjs::enable("adjustColorbinning")
  #               updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #             else{
  #               shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #               updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))}
  #           }
  #         }
  #
  #
  #         ###If Expression Assays##########################################################
  #       }else{Dfassay <- assay(vals$counts, input$AdvancedMethodSelect_Colorby)
  #       if(input$GeneSelect_Assays_Colorby %in% rownames(Dfassay)){
  #         Dfassay <- data.frame(Dfassay[which(rownames(Dfassay)== input$GeneSelect_Assays_Colorby),])
  #
  #         if(!is.numeric(Dfassay[,1])){
  #           updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #           shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #           shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #           updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #         }else if(is.integer(Dfassay[,1])
  #           &length(levels(as.factor(Dfassay[,1])))<=25
  #           &input$SelectColorType == 'Categorical'){
  #           updateCheckboxInput(session,"checkColorbinning","Perform Binning", value = FALSE)
  #           shinyjs::delay(5,shinyjs::disable("checkColorbinning"))
  #           shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #           updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))
  #
  #         }else if(is.integer(Dfassay[,1])
  #           &length(levels(as.factor(Dfassay[,1])))<=25
  #           &input$SelectColorType == 'Continuous'){
  #
  #           shinyjs::enable("checkColorbinning")
  #           if(input$checkColorbinning == TRUE){
  #             shinyjs::enable("adjustColorbinning")
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #           else{
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))}
  #
  #         }else{
  #
  #           shinyjs::enable("checkColorbinning")
  #           if(input$checkColorbinning == TRUE){
  #             shinyjs::enable("adjustColorbinning")
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("ggplot","Celda"))}
  #
  #           else{
  #             shinyjs::delay(5,shinyjs::disable("adjustColorbinning"))
  #             updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:", choices = c("RdYlBu",color_seqdiv))
  #           }
  #         }
  #       }
  #       }#Dfassay_end
  #     }#ifnot_end
  #   }
  # })###observe_end
  #


  #-+-+-+-+-+-cellviewer prepare step1: choose data. (next steps included)###########################################################
  cellviewer <- eventReactive(input$runCellViewer,{
    if(input$QuickAccess == ""){

    }else if(input$QuickAccess != "Custom"){
      ###QuickAccess for ReduceData
      xy <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$QuickAccess))
      colnames(xy) <- c("X_input","Y_input")
      xy <- cbind(xy,data.frame(colData(vals$counts)))

    }else{
      ###Custom
      #X_axis
      ##ReduceDim
      if(input$TypeSelect_Xaxis == "Reduced Dimensions"){
        Dfx <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Xaxis))
        Dfx2 <- Dfx[which(colnames(Dfx) == input$ColumnSelect_Xaxis)]
        colnames(Dfx2) <- c("X_input")
        rm(Dfx)##Assay
      }else if(input$TypeSelect_Xaxis == "Expression Assays"){
        Dfx <- assay(vals$counts, input$AdvancedMethodSelect_Xaxis)
        Dfx2 <- data.frame(Dfx[which(rownames(Dfx)== input$GeneSelect_Assays_Xaxis),])
        colnames(Dfx2) <- c("X_input")
        rm(Dfx)##Annotation
      }else if(input$TypeSelect_Xaxis == "Cell Annotation"){
        Dfx <- colData(vals$counts)
        Dfx2 <- data.frame(Dfx[which(colnames(Dfx)== input$AnnotationSelect_Xaxis)])
        colnames(Dfx2) <- c("X_input")
        rm(Dfx)
      }

      #Y_axis
      ##ReduceDIm
      if(input$TypeSelect_Yaxis == "Reduced Dimensions"){
        Dfy <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Yaxis))
        Dfy2 <- Dfy[which(colnames(Dfy) == input$ColumnSelect_Yaxis)]
        colnames(Dfy2) <- c("Y_input")
        rm(Dfy)##Assay
      }else if(input$TypeSelect_Yaxis == "Expression Assays"){
        Dfy <- assay(vals$counts, input$AdvancedMethodSelect_Yaxis)
        Dfy2 <- data.frame(Dfy[which(rownames(Dfy)== input$GeneSelect_Assays_Yaxis),])
        colnames(Dfy2) <- c("Y_input")
        rm(Dfy)##Annotation
      }else if(input$TypeSelect_Yaxis == "Cell Annotation"){
        Dfy <- colData(vals$counts)
        Dfy2 <- data.frame(Dfy[which(colnames(Dfy)== input$AnnotationSelect_Yaxis)])
        colnames(Dfy2) <- c("Y_input")
        rm(Dfy)
      }
      xy <- cbind(Dfx2,Dfy2)#BindXY
      xy <- cbind(xy,data.frame(colData(vals$counts)))#BindAnnotation
      rm(Dfx2)
      rm(Dfy2)
    }#ConditionalCustom_end

    #-+-+-+-+-+-cellviewer prepare2 : choose color#####################

    ####Cell Annotation if numeric, Categorical, check###

    if(input$TypeSelect_Colorby != 'Pick a Color'){
      ####Cell Annotation if numeric, Categorical, check###
      if(input$TypeSelect_Colorby == 'Cell Annotation'){
        if(!is.numeric(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])){
          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          # legendname <- paste0(input$AnnotationSelect_Colorby)
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
          &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25
          &input$SelectColorType == 'Categorical'){
          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          colnames(total_colors) <- c("Color")
          total_colors$Color <- as.factor(total_colors$Color)
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
          &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == FALSE){

          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
          &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)

        }else if(is.numeric(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)
        }else{

          total_colors <- colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)
        }

        ####Reduced Dimensions if numeric, Categorical, check###
      }else if(input$TypeSelect_Colorby == 'Reduced Dimensions'){
        Dfcolor <- data.frame(reducedDims(vals$counts)@listData[[input$ApproachSelect_Colorby]])
        Dfcolor <- Dfcolor[which(colnames(Dfcolor) == input$ColumnSelect_Colorby)]

        if(!is.numeric(Dfcolor[,1])){
          total_colors <- Dfcolor[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfcolor[,1])
          &length(levels(as.factor(Dfcolor[,1])))<=25
          &input$SelectColorType == 'Categorical'){
          total_colors <- Dfcolor[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          total_colors$Color <- as.factor(total_colors$Color)
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfcolor[,1])
          &length(levels(as.factor(Dfcolor[,1])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == FALSE){

          total_colors <- Dfcolor[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfcolor[,1])
          &length(levels(as.factor(Dfcolor[,1])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- Dfcolor[,1] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)

        }else if(is.numeric(Dfcolor[,1])
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- Dfcolor[,1] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)

        }else{

          total_colors <- Dfcolor[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)
        }

      }else{
        Dfassay <- assay(vals$counts, input$AdvancedMethodSelect_Colorby)
        Dfassay <- data.frame(Dfassay[which(rownames(Dfassay)== input$GeneSelect_Assays_Colorby),])

        if(!is.numeric(Dfassay[,1])){
          total_colors <- Dfassay[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfassay[,1])
          &length(levels(as.factor(Dfassay[,1])))<=25
          &input$SelectColorType == 'Categorical'){
          total_colors <- Dfassay[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          total_colors$Color <- as.factor(total_colors$Color)
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfassay[,1])
          &length(levels(as.factor(Dfassay[,1])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == FALSE){

          total_colors <- Dfassay[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)

        }else if(is.integer(Dfassay[,1])
          &length(levels(as.factor(Dfassay[,1])))<=25
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- Dfassay[,1] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)

        }else if(is.numeric(Dfassay[,1])
          &input$SelectColorType == 'Continuous'
          &input$checkColorbinning == TRUE){

          total_colors <- Dfassay[,1] %>% data.frame()
          color1 <- cut(total_colors[,1], breaks = seq(from = min(total_colors)-1,
            to = max(total_colors)+1,
            by = (max(total_colors)-min(total_colors)+1)/input$adjustColorbinning)) %>% data.frame()
          colnames(color1) <- c("Color")
          xy <- cbind(xy,color1)
          rm(color1)
          rm(total_colors)

        }else{

          total_colors <- Dfassay[,1] %>% data.frame()
          colnames(total_colors) <- c("Color")
          xy <- cbind(xy,total_colors)
          rm(total_colors)
        }

      }

    }#ifnotUniform_end

    #-+-+-+-+-+-cellviewer prepare3 : prepare Axis Label Name#####################
    ###Xaxis label name
    if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustxlab == ""){
      xname = paste0(input$QuickAccess, 1)
    }else if(input$QuickAccess != "Custom" & input$QuickAccess != ""& input$adjustxlab != ""){
      xname = input$adjustxlab
    }else if(input$TypeSelect_Xaxis == 'Reduced Dimensions'){
      xname = paste0(input$ApproachSelect_Xaxis,substr(input$ColumnSelect_Xaxis,2,2))
    }else if(input$TypeSelect_Xaxis == 'Expression Assays'){
      xname = paste0(input$GeneSelect_Assays_Xaxis)
    }else{
      xname = paste0(input$AnnotationSelect_Xaxis)
    }

    ###Yaxis label name
    if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustylab == ""){
      yname = paste0(input$QuickAccess, 2)
    }else if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustylab != ""){
      yname = input$adjustylab
    }else if(input$TypeSelect_Yaxis == 'Reduced Dimensions'){
      yname = paste0(input$ApproachSelect_Yaxis,substr(input$ColumnSelect_Yaxis,2,2))
    }else if(input$TypeSelect_Yaxis == 'Expression Assays'){
      yname = paste0(input$GeneSelect_Assays_Yaxis)

    }else{
      yname = paste0(input$AnnotationSelect_Yaxis)
    }

    ###Yaxis label name
    if(input$TypeSelect_Colorby != 'Pick a Color'){
      if(input$TypeSelect_Colorby == 'Reduced Dimensions'){
        legendname = paste0(input$ApproachSelect_Colorby,substr(input$ColumnSelect_Colorby,2,2))

      }else if(input$TypeSelect_Colorby == 'Expression Assays'){
        legendname = paste0(input$GeneSelect_Assays_Colorby)

      }else{
        legendname = paste0(input$AnnotationSelect_Colorby)
      }
    }

    #-+-+-+-+-+-cellviewer prepare4 : choose group by and create plotly function###################

    if (input$adjustgroupby == "None"){
      #if uniform
      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() + xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }
      #if not uniform
      else{
        #ggplot#none
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() + xlab(xname) + ylab(paste0("\n",yname)) +  labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }
        #ggplotly#none
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600) }
      #else_end

    }#if_none_end

    ###Integer,level>25
    else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
      & length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))>25){
      #data manage#Integer,level>25
      total_features <- colData(vals$counts)@listData[[input$adjustgroupby]]
      c1 <- cut(total_features, breaks = seq(from = min(total_features)-1,
        to = max(total_features)+1,
        by = (max(total_features)-min(total_features)+1)/input$adjustbinning)) %>%
        data.frame()
      colnames(c1) <- c("groupby")
      c1$groupby <- as.factor(c1$groupby)
      xy <- cbind(xy,c1)
      rm(c1)

      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }
      else{
        #ggplot#Integer,level>25
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname)) + labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }
        #ggplotly#Integer,level>25
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600)
      }
    }

    ###Integer,level<25,continuous
    else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
      &length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))<=25
      &input$SelectValueType == "Continuous"){
      #data manage#Integer,level<25,Continuous
      total_features <- colData(vals$counts)@listData[[input$adjustgroupby]]
      c1 <- cut(total_features, breaks = seq(from = min(total_features)-1,
        to = max(total_features)+1,
        by = (max(total_features)-min(total_features)+1)/input$adjustbinning)) %>%
        data.frame()
      colnames(c1) <- c("groupby")
      c1$groupby <- as.factor(c1$groupby)
      xy <- cbind(xy,c1)
      rm(c1)

      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }#ifUniform_end
      else{
        #ggplot#Integer,level<25,Continous
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname)) + labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }

        #ggplotly#Integer,level<25,Continous
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600)
      }#notUniform_End
    }#condition_End

    ###Integer,level<25,Categorical
    else if(is.integer(colData(vals$counts)@listData[[input$adjustgroupby]])
      &length(levels(as.factor(colData(vals$counts)@listData[[input$adjustgroupby]])))<=25
      &input$SelectValueType == "Categorical"){
      #data manage#Integer,level<25,Categorical
      c1 <- colData(vals$counts)@listData[[input$adjustgroupby]] %>% data.frame()
      colnames(c1) <- c("groupby")
      c1$groupby <- as.factor(c1$groupby)
      xy <- cbind(xy,c1)
      rm(c1)

      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }#uniform_end
      else{
        #ggplot#Integer,level<25,Categorical
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname)) + labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }
        #ggplotly#Integer,level<25,Categorical
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600)
      }#notuniform_End
    }#condition_End

    ###Numeric,noninteger
    else if (is.numeric(colData(vals$counts)@listData[[input$adjustgroupby]])){
      #data manage#Numeric,noninteger
      total_features <- colData(vals$counts)@listData[[input$adjustgroupby]]
      c1 <- cut(total_features,
        breaks = seq(from = min(total_features)-1, to = max(total_features)+1,
          by = (max(total_features)-min(total_features)+1)/input$adjustbinning)) %>% data.frame()
      colnames(c1) <- c("groupby")
      c1$groupby <- as.factor(c1$groupby)
      xy <- cbind(xy,c1)
      rm(c1)

      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }#ifUniform_end
      else{
        #ggplot2#Numeric,noninteger
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname)) + labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }
        #ggplotly2#Numeric,noninteger
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600)
      }#notUniform_end
    }#condition_end

    ###else,Categorical
    else{
      #data manage#cate
      c1 <- colData(vals$counts)@listData[[input$adjustgroupby]] %>% data.frame()
      colnames(c1) <- c("groupby")
      c1$groupby <- as.factor(c1$groupby)
      xy <- cbind(xy,c1)
      rm(c1)
      if(input$TypeSelect_Colorby == 'Pick a Color'){
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input") +
          geom_point(color = input$Col, size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname))
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input"), height = 600)
      }#ifUniform_end
      else{
        #ggplot3#
        a <- ggplot(data = xy) +
          aes_string(x= "X_input", y= "Y_input", color = "Color") +
          geom_point(size = input$adjustsize, alpha = input$adjustalpha) +
          theme_classic() +
          theme(legend.title = element_blank(),
            strip.background = element_blank()) +
          facet_wrap(~groupby) +
          xlab(xname) + ylab(paste0("\n",yname)) + labs(color= legendname)

        if(!is.numeric(xy$Color)){
          if(input$adjustbrewer == 'Celda'){
            a = a + scale_color_manual(values = celda::distinctColors(length(levels(xy$Color)))) + theme(legend.text=element_text(size=12))}
          else{a = a + theme(legend.text=element_text(size=12))}
        }else{
          a = a + scale_color_distiller(palette = input$adjustbrewer)
        }
        #ggplotly3#
        if (input$adjusttitle != ""){
          a <- a + ggtitle(input$adjusttitle)
        }
        ggplotly(a, tooltip = c("X_input", "Y_input", "Color"), height = 600)

      }#notUniform_end
    }#condition_end

  })#Cellviewer_end
  output$scatter <- renderPlotly({cellviewer()})
  #
  #
  #-+-+-+-+-+-cellviewer prepare done: plot#####################
  ###plotly_after_reactive


  #-----------------------------------------------------------------------------
  # Page 4: Batch Correction ####
  #-----------------------------------------------------------------------------

  observeEvent(input$toggleNormalization, {
    if (vals$showAssayDetails == FALSE) {
      vals$showAssayDetails <- TRUE
      shinyjs::show(id="normalization", anim=TRUE, animType="slide", time=0.2)
      updateActionButton(session, "toggleAssayDetails", icon=icon("caret-up", lib="font-awesome"))
    } else {
      vals$showAssayDetails <- FALSE
      shinyjs::hide(id="normalization", anim=TRUE, animType="slide", time=0.2)
      updateActionButton(session, "toggleAssayDetails", icon=icon("caret-down", lib="font-awesome"))
    }
  })

  output$batchCheckResUI <- renderUI({
    if(input$batchCheckResType == 1){
      selectInput("batchCheckCorrAssay", "Corrected Assay", c("", vals$batchResAssay))
    } else {
      selectInput("batchCheckCorrReddim", "Corrected Reduced Dimension", c("", vals$batchResReddim))
    }
  })

  output$batchOriVar <- renderPlot({
    if (!is.null(vals$counts) &
        input$batchCheckVar != "None"){
      plotSCEBatchFeatureMean(inSCE = vals$counts,
        useAssay = input$batchCheckOrigAssay,
        batch = input$batchCheckVar)
    }
  })

  output$batchOriPCA <- renderPlot({
    if(!is.null(vals$counts) &
        input$batchCheckVar != "None" &
        input$batchCheckVar != input$batchCheckCond){
      if(input$batchCheckCond == "None"){
        shapeBy <- NULL
      } else {
        shapeBy <- input$batchCheckCond
      }
      pcaName <- paste0(input$batchCheckOrigAssay, "_PCA")
      if(!"PCA" %in% names(SingleCellExperiment::reducedDims(vals$counts))){
        vals$counts <- getPCA(vals$counts, useAssay = input$batchCheckOrigAssay,
          reducedDimName = pcaName)
        updateReddimInputs()
      }
      plotSCEDimReduceColData(vals$counts, colorBy = input$batchCheckVar,
        shape = shapeBy, reducedDimName = pcaName,
        title = paste0("Original ", input$batchCheckOrigAssay, " PCA"))
    }
  })

  output$batchCorrVar <- renderPlot({
    if (!is.null(vals$counts) &
        input$batchCheckVar != "None"){
      if(input$batchCheckResType == 1 &
          length(vals$batchResAssay) > 0 &
          !is.null(input$batchCheckCorrAssay)){
        if(input$batchCheckCorrAssay != ""){
          plotSCEBatchFeatureMean(inSCE = vals$counts,
            useAssay = input$batchCheckCorrAssay,
            batch = input$batchCheckVar)
        }
      } else if(input$batchCheckResType == 2 &
          length(vals$batchResReddim) > 0 &
          !is.null(input$batchCheckCorrReddim)){
        if(input$batchCheckCorrReddim != ""){
          plotSCEBatchFeatureMean(inSCE = vals$counts,
            useReddim = input$batchCheckCorrReddim,
            batch = input$batchCheckVar)
        }
      }
    }
  })

  output$batchCorrReddim <- renderPlot({
    if (!is.null(vals$counts) &
        input$batchCheckVar != "None" &
        input$batchCheckVar != input$batchCheckCond){
      if(input$batchCheckCond == "None"){
        shapeBy = NULL
      } else {
        shapeBy = input$batchCheckCond
      }
      if(input$batchCheckResType == 2 &
          length(vals$batchResReddim) > 0 &
          !is.null(input$batchCheckCorrReddim)){
        if(input$batchCheckCorrReddim != ""){
          plotSCEDimReduceColData(vals$counts, colorBy = input$batchCheckVar,
            shape = shapeBy, reducedDimName = input$batchCheckCorrReddim,
            conditionClass = "character",
            title = paste0(input$batchCheckCorrReddim, " corrected"))
        }
      } else if (input$batchCheckResType == 1 &
          length(vals$batchResAssay) > 0 &
          !is.null(input$batchCheckCorrAssay)){
        if(input$batchCheckCorrAssay != ""){
          pcaName <- paste0(input$batchCheckCorrAssay, "_PCA")
          vals$counts <- getPCA(vals$counts, useAssay = input$batchCheckCorrAssay,
            reducedDimName = pcaName)
          updateReddimInputs()
          plotSCEDimReduceColData(vals$counts, colorBy = input$batchCheckVar,
            shape = shapeBy, reducedDimName = pcaName,
            title = paste0(input$batchCheckCorrAssay, " corrected"))
        }
      }
    }
  })

  observeEvent(input$BBKNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("BBKNNRun", {
        saveassayname <- gsub(" ", "_", input$BBKNNSaveReddim)
        vals$counts <- runBBKNN(vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          reducedDimName = saveassayname,
          nComponents = input$BBKNNNComp)
        shinyalert::shinyalert('Success!', 'BBKNN completed.',
          type = 'success')
        vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
        updateReddimInputs()
        vals$batchCorrStatus <- "BBKNN Complete"
      }
      )
    }
  })

  output$selectCombatRefBatchUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$combatRef){
        selectInput("combatRefBatch", "Choose Reference Batch:",
          unique(sort(colData(vals$counts)[, input$batchCorrVar])))
      }
    }
  })

  observeEvent(input$combatRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("combatRun", {
        #check for zeros
        if (any(rowSums(assay(vals$counts, input$batchCorrAssay)) == 0)){
          shinyalert::shinyalert("Error!", "Rows with a sum of zero found. Filter data to continue.", type = "error")
        }else {
          saveassayname <- gsub(" ", "_", input$combatSaveAssay)
          if (input$combatRef){
            assay(vals$counts, saveassayname) <-
              ComBatSCE(inSCE = vals$counts, batch = input$batchCorrVar,
                useAssay = input$batchCorrAssay,
                par.prior = input$combatParametric,
                covariates = input$combatCond,
                mean.only = input$combatMeanOnly,
                ref.batch = input$combatRefBatch)
          } else {
            assay(vals$counts, saveassayname) <-
              ComBatSCE(inSCE = vals$counts, batch = input$batchCorrVar,
                useAssay = input$batchCorrAssay,
                par.prior = input$combatParametric,
                covariates = input$combatCond,
                mean.only = input$combatMeanOnly)
          }
          vals$batchResAssay <- c(vals$batchResAssay, saveassayname)
          updateAssayInputs()
          shinyalert::shinyalert('Success!', 'ComBat completed.', type = 'success')
          vals$batchCorrStatus <- "ComBat Complete"

        }
      })
    }
  })

  observeEvent(input$FastMNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("FastMNNRun", {
        saveassayname <- gsub(" ", "_", input$FastMNNSaveReddim)
        if(isTRUE(input$FastMNNPcInput)){
          fmnnAssay <- input$FastMNNReddim
        } else {
          fmnnAssay <- input$batchCorrAssay
        }
        vals$counts <- runFastMNN(vals$counts,
          useAssay = fmnnAssay,
          batch = input$batchCorrVar,
          reducedDimName = saveassayname,
          pcInput = input$FastMNNPcInput
        )
        shinyalert::shinyalert('Success!', 'FastMNN completed.',
          type = 'success')
        vals$batchCorrStatus <- "FastMNN Complete"
        vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
        updateReddimInputs()
      }
      )
    }
  })

  observeEvent(input$HarmonyRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("HarmonyRun", {
        saveassayname <- gsub(" ", "_", input$HarmonySaveReddim)
        if(isTRUE(input$HarmonyPcInput)){
          useAssay <- input$HarmonyReddim
        } else {
          useAssay <- input$batchCorrAssay
        }
        if(is.na(as.numeric(input$HarmonyTheta))){
          vals$batchCorrStatus <- ""
          stop("Theta value must be numeric.")
        } else {
          theta <- as.numeric(input$HarmonyTheta)
        }
        vals$counts <- runHarmony(vals$counts, useAssay = useAssay,
          pcInput = input$HarmonyPcInput,
          batch = input$batchCorrVar,
          reducedDimName = saveassayname,
          nComponents = input$HarmonyNComp,
          theta = theta, nIter = input$HarmonyNIter)
        shinyalert::shinyalert('Success!', 'Harmony completed.',
          type = 'success')
        vals$batchCorrStatus <- "Harmony Complete"
        vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
        updateReddimInputs()
      }
      )
    }
  })

  observeEvent(input$limmaRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("limmaRun", {
        saveassayname <- gsub(" ", "_", input$limmaSaveAssay)
        vals$counts <- runLimmaBC(vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          assayName = saveassayname)
        shinyalert::shinyalert('Success!', 'Limma completed.',
          type = 'success')
        vals$batchCorrStatus <- "Limma Complete"
        vals$batchResAssay <- c(vals$batchResAssay, saveassayname)
        updateAssayInputs()
      }
      )
    }
  })

  observeEvent(input$ligerRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("ligerRun", {
        #check for zeros
        if (any(rowSums(assay(vals$counts, input$batchCorrAssay)) == 0)){
          shinyalert::shinyalert("Error!", "Rows with a sum of zero found. Filter data to continue.", type = "error")
        } else {
          saveassayname <- gsub(" ", "_", input$ligerSaveReddim)
          vals$counts <-
            runLIGER(inSCE = vals$counts,
              useAssay = input$batchCorrAssay,
              batch = input$batchCorrVar,
              reducedDimName = saveassayname,
              nComponents = input$ligerNComp,
              lambda = input$ligerLambda,
              resolution = input$ligerResolution)
          shinyalert::shinyalert('Success!', 'LIGER completed.',
            type = 'success')
          vals$batchCorrStatus <- "LIGER Complete"
          vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
          updateReddimInputs()
        }
      })
    }
  })

  observeEvent(input$MNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("MNNRun", {
        saveassayname <- gsub(" ", "_", input$MNNSaveAssay)
        if(is.na(as.numeric(input$MNNSigma))){
          vals$batchCorrStatus <- ""
          stop("Sigma value must be numeric.")
        } else {
          sigma <- as.numeric(input$MNNSigma)
        }
        vals$counts <- runMNNCorrect(vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          k = input$MNNK, sigma = sigma,
          assayName = saveassayname)
        shinyalert::shinyalert('Success!', 'MNN completed.',
          type = 'success')
        vals$batchCorrStatus <- "MNN Complete"
        vals$batchResAssay <- c(vals$batchResAssay, saveassayname)
        updateAssayInputs()
      }
      )
    }
  })

  output$scMergeNBatch <- renderUI({
    if(!is.null(vals$counts) &&
        !is.null(input$batchCorrVar) &&
        !input$batchCorrVar == 'None'){
      nBatch <- length(unique(SummarizedExperiment::colData(vals$counts)[[input$batchCorrVar]]))
      span(paste0("Please input ", nBatch, " integer(s), separated by ','."))
    }
  })

  observeEvent(input$scMergeRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("scMergeRun", {
        saveassayname <- gsub(" ", "_", input$scMergeSaveAssay)
        if(input$scMergeSEGOpt == 1){
          seg <- NULL
        } else if(input$scMergeSEGOpt == 2){
          data("SEG")
          seg <- SEG[[input$scMergeSEGSpecies]]
        } else {
          seg <- str_trim(scan(text = input$scMergeSEGCustom,
            sep='\n', what = 'character'))
        }
        if(isTRUE(input$scMergeAutoKmk)){
          kmk <- NULL
        } else {
          kmk <- scan(text = input$scMergeUserKmk, sep=',')
        }
        vals$counts <- runSCMerge(inSCE = vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          cellType = input$scMergeCT,
          seg = seg, kmeansK = kmk,
          assayName = saveassayname
        )
        shinyalert::shinyalert('Success!', 'scMerge completed.',
          type = 'success')
        vals$batchCorrStatus <- "scMerge Complete"
        vals$batchResAssay <- c(vals$batchResAssay, saveassayname)
        updateAssayInputs()
      })
    }
  })

  output$Srt3IntNAnchUI <- renderUI({
    if(!is.null(vals$counts)){
      ngene <- nrow(vals$counts)
      numericInput('Srt3IntNAnch', "Number of anchors",
        value = ngene, min = 30, max = ngene, step = 1)
    }
  })

  observeEvent(input$Srt3IntRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("Srt3IntRun", {
        saveassayname <- gsub(" ", "_", input$Srt3IntSaveAssay)
        vals$counts <- runSeurat3Integration(vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          assayName = saveassayname,
          nAnchors = input$Srt3IntNAnch
        )
        # According to input nAnchor, corrected matrix can be an "assay" or a
        # "reducedDim". Here I temprorily make a condition basing on how Seurat
        # performs empirically.
        if(input$Srt3IntNAnch == nrow(vals$counts)){
          # Usually in this condition, seurat returns a full-sized assay
          vals$batchResAssay <- c(vals$batchResAssay, saveassayname)
          updateAssayInputs()
        } else if(input$Srt3IntNAnch == nrow(vals$counts)){
          # Under this condition, seurat usually returns a reduced matrix of
          # <= nAnchor dimensions.
          vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
          updateReddimInputs()
        }
        shinyalert::shinyalert('Success!', 'Seurat3 Integration completed.',
          type = 'success')
        vals$batchCorrStatus <- "Seurat3 Integration Complete"
      }
      )
    }
  })

  output$zinbwaveNHvgUI <- renderUI({
    if(!is.null(vals$counts)){
      ngenes <- nrow(vals$counts)
      zwdefault <- min(ngenes, 1000)
      numericInput('zinbwaveNHVG', 'Number of highly variable genes to use:',
        value = zwdefault, max = ngenes)
    }
  })

  output$zinbwaveEpsUI <- renderUI({
    if(!is.null(vals$counts)){
      ngenes <- nrow(vals$counts)
      zwdefault <- min(ngenes, 1000)
      numericInput('zinbwaveEps', 'Epsilon value:',
        value = zwdefault, max = ngenes)
    }
  })

  observeEvent(input$zinbwaveRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("zinbwaveRun", {
        saveassayname <- gsub(" ", "_", input$limmaSaveAssay)
        vals$counts <- runZINBWaVE(vals$counts,
          useAssay = input$batchCorrAssay,
          batch = input$batchCorrVar,
          reducedDimName = saveassayname,
          epsilon = input$zinbwaveEps,
          nHVG = input$zinbwaveNHVG,
          nIter = input$zinbwaveNIter,
          nComponents = input$zinbwaveNComp
        )
        shinyalert::shinyalert('Success!', 'ZINBWaVE completed.',
          type = 'success')
        vals$batchCorrStatus <- "ZINBWaVE Complete"
        vals$batchResReddim <- c(vals$batchResReddim, saveassayname)
        updateReddimInputs()
      }
      )
    }
  })

  output$batchCorrStatus <- renderUI({
    span(vals$batchCorrStatus, style = "color:green;margin-top:30px;")
  })

  #-----------------------------------------------------------------------------
  # Page 4.1: Feature Selection
  #-----------------------------------------------------------------------------
    observeEvent(input$findHvgButtonFS, {
        if (!is.null(vals$counts)) {
            if (input$hvgMethodFS == "vst"
                || input$hvgMethodFS == "mean.var.plot"
                || input$hvgMethodFS == "dispersion") {
                withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
                    #vals$counts <- seuratFindHVG(vals$counts, useAssay = input$assaySelectFS, seuratWorkflow$geneNamesSeurat, input$hvgMethodFS, as.numeric(input$hvgNoFeaturesFS))
                    vals$counts <- seuratFindHVG(inSCE = vals$counts,
                                                 useAssay = input$assaySelectFS,
                                                 hvgMethod = input$hvgMethodFS,
                                                 hvgNumber = as.numeric(input$hvgNoFeaturesFS))
                    if (input$hvgMethodFS == "vst") {
                        vals$vfplot <- ggplot() + geom_point(aes(x = log(rowData(vals$counts)$seurat_variableFeatures_vst_mean), y = rowData(vals$counts)$seurat_variableFeatures_vst_varianceStandardized)) + geom_point(aes(x = log(subset(rowData(vals$counts)$seurat_variableFeatures_vst_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS)))), y = subset(rowData(vals$counts)$seurat_variableFeatures_vst_varianceStandardized, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS)))), colour = "red") + geom_label(aes(x = log(subset(rowData(vals$counts)$seurat_variableFeatures_vst_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS)))), y = subset(rowData(vals$counts)$seurat_variableFeatures_vst_varianceStandardized, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), label = subset(rownames(vals$counts), rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS)))), colour = "red", size = 2) + labs(x = "Mean", y = "Standardized Variance")
                    }
                    else if (input$hvgMethodFS == "mean.var.plot") {
                        vals$vfplot <- ggplot() + geom_point(aes(x = rowData(vals$counts)$seurat_variableFeatures_mvp_mean, y = rowData(vals$counts)$seurat_variableFeatures_mvp_dispersionScaled)) + geom_point(aes(x = subset(rowData(vals$counts)$seurat_variableFeatures_mvp_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS))), y = subset(rowData(vals$counts)$seurat_variableFeatures_mvp_dispersionScaled, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS)))), colour = "red") + geom_label(aes(x = subset(rowData(vals$counts)$seurat_variableFeatures_mvp_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), y = subset(rowData(vals$counts)$seurat_variableFeatures_mvp_dispersionScaled, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), label = subset(rownames(vals$counts), rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS)))), colour = "red", size = 2) + labs(x = "Mean", y = "Dispersion")
                    }
                    else if (input$hvgMethodFS == "dispersion") {
                        vals$vfplot <- ggplot() + geom_point(aes(x = rowData(vals$counts)$seurat_variableFeatures_dispersion_mean, y = rowData(vals$counts)$seurat_variableFeatures_dispersion_dispersionScaled)) + geom_point(aes(x = subset(rowData(vals$counts)$seurat_variableFeatures_dispersion_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS))), y = subset(rowData(vals$counts)$seurat_variableFeatures_dispersion_dispersionScaled, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS)))), colour = "red") + geom_label(aes(x = subset(rowData(vals$counts)$seurat_variableFeatures_dispersion_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), y = subset(rowData(vals$counts)$seurat_variableFeatures_dispersion_dispersionScaled, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), label = subset(rownames(vals$counts), rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS)))), colour = "red", size = 2) + labs(x = "Mean", y = "Dispersion")
                    }
                })
                showNotification("Find HVG Complete")
            }
            else if (input$hvgMethodFS == "modelGeneVar") {
                withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
                    vals$counts <- scran_modelGeneVar(inSCE = vals$counts, assayName = input$assaySelectFS)
                    vals$vfplot <- ggplot() + geom_point(aes(x = rowData(vals$counts)$scran_modelGeneVar_mean, y = rowData(vals$counts)$scran_modelGeneVar_totalVariance)) + geom_point(aes(x = subset(rowData(vals$counts)$scran_modelGeneVar_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS))), y = subset(rowData(vals$counts)$scran_modelGeneVar_totalVariance, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesFS)))), colour = "red") + geom_label(aes(x = subset(rowData(vals$counts)$scran_modelGeneVar_mean, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), y = subset(rowData(vals$counts)$scran_modelGeneVar_totalVariance, rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))), label = subset(rownames(vals$counts), rownames(vals$counts) %in% getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS)))), colour = "red", size = 2) + labs(x = "Average Expression", y = "Variance")
                })
                showNotification("Scran modelGeneVar processing complete!")
            }
        }
        else if (is.null(vals$counts)) {
            showNotification("Please input dataset (rds file) before computing highly variable genes!", type = "error")
        }
        else {
            showNotification("An error occurred while computing highly variable genes!", type = "error")
        }
    })

    output$plotFS <- renderPlot({
        if (!is.null(vals$hvgPlotFS)) {
            vals$hvgPlotFS
        }
    })

    observe({
        if (!is.null(vals$vfplot)) {
            if (!is.na(as.numeric(input$hvgNoFeaturesViewFS))) {
                vals$hvgPlotFS <- vals$vfplot
            }
        }
    })

    output$hvgOutputFS <- renderText({
    if (!is.null(vals$counts)) {
      if (!is.null(vals$vfplot)) {
        getTopHVG(inSCE = vals$counts, method = input$hvgMethodFS, n = as.numeric(input$hvgNoFeaturesViewFS))
      }
    }
 })

  #-----------------------------------------------------------------------------
  # Page 5.1: Differential Expression
  #-----------------------------------------------------------------------------
  shinyjs::onclick("Diffex_hideAllSections", allSections(
    "hide", c(paste("de", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("Diffex_showAllSections", allSections(
    "show", c(paste("de", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("diffex1",
                   shinyjs::toggle(id = "de1",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex2",
                   shinyjs::toggle(id = "de2",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex3",
                   shinyjs::toggle(id = "de3",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex4",
                   shinyjs::toggle(id = "de4",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex5",
                   shinyjs::toggle(id = "de5",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex6",
                   shinyjs::toggle(id = "de6",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex7",
                   shinyjs::toggle(id = "de7",
                                   anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "diffex1", class = "btn-block")
  shinyjs::addClass(id = "diffex2", class = "btn-block")
  shinyjs::addClass(id = "diffex3", class = "btn-block")
  shinyjs::addClass(id = "diffex4", class = "btn-block")
  shinyjs::addClass(id = "diffex5", class = "btn-block")
  shinyjs::addClass(id = "diffex6", class = "btn-block")
  shinyjs::addClass(id = "diffex7", class = "btn-block")

  output$selectDiffexConditionUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$selectDiffex == "ANOVA") {
        tagList(
          selectInput("selectDiffexCondition", "Select Condition(s):",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      } else {
        tagList(
          selectInput("selectDiffexCondition",
                      "Select Condition:",
                      colnames(colData(vals$counts))),
          selectInput("selectDiffexCovariates",
                      "Select Additional Covariates:",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      }
    }
  })

  #For conditions with more than two factors, select the factor of interest
  output$selectDiffexConditionLevelUI <- renderUI({
    req(vals$counts)
    if (length(colnames(colData(vals$counts))) > 0){
      if (length(unique(colData(vals$counts)[, input$selectDiffexCondition])) > 2 & input$selectDiffex == "DESeq2"){
        tagList(
          radioButtons("selectDiffexConditionMethod", "Select Analysis Method:",
                       choiceNames = c("Biomarker (1 vs all)", "Factor of Interest vs. Control Factor",
                                       "Entire Factor (Full/Reduced)"),
                       choiceValues = c("biomarker", "contrast", "fullreduced")
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod != 'fullreduced'",
            selectInput("selectDiffexConditionOfInterest",
                        "Select Factor of Interest",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod == 'contrast' && input.selectDiffex == 'DESeq2'",
            selectInput("selectDiffexControlCondition",
                        "Select Control Factor",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
          )
        )
      } else if (length(unique(colData(vals$counts)[, input$selectDiffexCondition])) > 2 & input$selectDiffex == "limma") {
        tagList(
          radioButtons("selectDiffexConditionMethod", "Select Analysis Method:",
                       choiceNames = c("Biomarker (1 vs all)", "Factor of Interest",
                                       "Entire Factor"),
                       choiceValues = c("biomarker", "coef", "allcoef")
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod != 'allcoef'",
            selectInput("selectDiffexConditionOfInterest",
                        "Select Factor of Interest",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
          )
        )
      } else if (input$selectDiffex == "ANOVA") {
        tagList(
          selectInput("anovaCovariates", "Select Additional Covariates:",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      }
    }
  })

  #Run differential expression
  observeEvent(input$runDiffex, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("runDiffex", {
        vals$diffexheatmapplot <- NULL
        #run diffex to get gene list and pvalues
        if (input$selectDiffex == "ANOVA"){
          useCovariates <- input$anovaCovariates
        } else {
          useCovariates <- input$selectDiffexCovariates
        }
        vals$diffexgenelist <- scDiffEx(inSCE = vals$counts,
                                        useAssay = input$diffexAssay,
                                        condition = input$selectDiffexCondition,
                                        covariates = useCovariates,
                                        significance = input$selectPval,
                                        ntop = nrow(vals$counts),
                                        usesig = FALSE,
                                        diffexmethod = input$selectDiffex,
                                        levelofinterest = input$selectDiffexConditionOfInterest,
                                        analysisType = input$selectDiffexConditionMethod,
                                        controlLevel = input$selectDiffexControlCondition,
                                        adjust = input$selectCorrection)
      })
    }
  })

  output$colorBarConditionUI <- renderUI({
    if (is.null(vals$counts)){
      selectInput("colorBarCondition", "Select Condition", NULL)
    } else {
      selectInput("colorBarCondition", "Select Condition",
                  colnames(colData(vals$counts)), multiple = TRUE)
    }
  })

  annotationColors <- reactiveValues(cols = list())

  output$heatmapSampleAnnotations <- renderUI({
    if (!is.null(input$colorBarCondition)) {
      if (!is.null(vals$counts) & length(input$colorBarCondition) > 0){
        if (all(input$colorBarCondition %in% colnames(colData(vals$counts)))) {
          h <- input$colorBarCondition
          L <- lapply(seq_along(h), function(i) colourGroupInput(paste0("colorGroup", i)))
          annotationColors$cols <- lapply(
            seq_along(h),
            function(i) {
              callModule(colourGroup, paste0("colorGroup", i), heading = h[i],
                         options = unique(unlist(colData(vals$counts)[, h[i]])))
            }
          )
          return(L)
        }
      }
    }
  })

  output$diffexNgenes <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Max genes: "), nrow(vals$diffexgenelist), sep = ""))
  })

  output$logFCDiffexRange <- renderUI({
    req(vals$diffexgenelist)
    if (input$selectDiffex != 'ANOVA') {
      logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
      if (length(logFCIndex) == 0) {
        #for DESeq2 with more than 1 covariate, choose the first column
        logFCIndex <- 1
      }
      minlogFC <- paste(em("Min logFC : "), round(min(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      maxlogFC <- paste(em("Max logFC : "), round(max(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      HTML(paste(minlogFC, maxlogFC, sep = '<br/>'))
    }
  })

  #Plot the differential expression results
  observeEvent(input$runPlotDiffex, {
    req(vals$diffexgenelist)
    withBusyIndicatorServer("runPlotDiffex", {
      tryCatch ({
        #logFC or abs(logFC)
        if (input$applyAbslogFCDiffex == TRUE) {
          absLogFCDiffex <- abs(input$selectlogFCDiffex)
        } else {
          absLogFCDiffex <- input$selectlogFCDiffex
        }
        #for convenience, index logFC and p-val columns for all the methods
        pvalIndex <- which(grepl("*padj*", colnames(vals$diffexgenelist)))
        logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
        if (input$selectNGenes > nrow(vals$diffexgenelist)) {
          stop("Max value exceeded for Input.")
        }
        #p-Val and logFC cutoff
        if (input$applyCutoff == TRUE & input$applylogFCCutoff == TRUE) {
          if (input$selectDiffex == 'ANOVA') {
            stop("logFC is not applicable for ANOVA")
          } else {
            if (min(na.omit(vals$diffexgenelist[, pvalIndex])) > input$selectPval) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least p-value in the results is greater than the selected p-val range")
            } else if (min(na.omit(vals$diffexgenelist[, logFCIndex])) > absLogFCDiffex) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least logFC in the results is greater than the selected logFC range")
            } else {
              diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, pvalIndex] <= input$selectPval &
                                                         vals$diffexgenelist[, logFCIndex] <= absLogFCDiffex), ]
            }
          }
        }
        #p-Val cutoff
        else if (input$applyCutoff == TRUE) {
          if (min(na.omit(vals$diffexgenelist[, pvalIndex])) > input$selectPval) {
            diffexFilterRes <- vals$diffexgenelist
            stop("the min/least p-value in the results is greater than the selected p-val range")
          } else {
            diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, pvalIndex] <= input$selectPval), ]
          }
        }
        #logFC cutoff
        else if (input$applylogFCCutoff == TRUE) {
          if (input$selectDiffex == 'ANOVA') {
            stop("logFC is not applicable for ANOVA")
          } else  {
            if (min(na.omit(vals$diffexgenelist[, logFCIndex])) > absLogFCDiffex) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least logFC in the results is greater than the selected logFC range")
            } else {
              diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, logFCIndex] <= absLogFCDiffex), ]
            }
          }
        } else {
          diffexFilterRes <- vals$diffexgenelist
        }
        if (is.null(diffexFilterRes)){
          diffexFilterRes <- vals$diffexgenelist
        }
        rowLengthFiltered <- nrow(diffexFilterRes)
        if (rowLengthFiltered == 0) {
          stop("You've got 0 genes after filtering.. adjust your filters accordingly")
        }
        if (rowLengthFiltered < input$selectNGenes) {
          diffexFilterRes <- diffexFilterRes[seq_len(rowLengthFiltered), ]
        } else {
          diffexFilterRes <- diffexFilterRes[seq_len(input$selectNGenes), ]
        }
        #run plotDiffex
        if (!is.null(vals$diffexgenelist)){
          if (input$displayHeatmapColorBar){
            if (is.null(input$colorBarCondition)){
              colors <- NULL
            } else {
              colors <- lapply(annotationColors$cols, function(col) col())
              names(colors) <- input$colorBarCondition
              if (is.null(colors[[length(colors)]][[1]])){
                colors <- NULL
              }
            }
          } else {
            colors <- NULL
          }
          vals$diffexheatmapplot <- plotDiffEx(inSCE = vals$counts,
                                               useAssay = input$diffexAssay,
                                               condition = input$colorBarCondition,
                                               geneList = rownames(diffexFilterRes),
                                               clusterRow = input$clusterRows,
                                               clusterCol = input$clusterColumns,
                                               displayRowLabels = input$displayHeatmapRowLabels,
                                               displayColumnLabels = input$displayHeatmapColumnLabels,
                                               displayRowDendrograms = input$displayHeatmapRowDendrograms,
                                               displayColumnDendrograms = input$displayHeatmapColumnDendrograms,
                                               annotationColors = colors,
                                               scaleExpression = input$applyScaleDiffex,
                                               columnTitle = input$heatmapColumnsTitle)
        }
      }, error = function(e){
        shinyalert::shinyalert("Error!", e$message, type = "error")
      })
    })
  })

  output$diffPlot <- renderPlot({
    req(vals$diffexheatmapplot)
    ComplexHeatmap::draw(vals$diffexheatmapplot)
  }, height = 600)

  #Create the differential expression results table
  output$diffextable <- DT::renderDataTable({
    if (!is.null(vals$diffexgenelist)){
      temptable <- cbind(rownames(vals$diffexgenelist), data.frame(vals$diffexgenelist))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, rownames = FALSE)

  #disable downloadGeneList button if the result is not null
  isDiffExResult <- reactive(is.null(vals$diffexgenelist))
  observe({
    if (isDiffExResult()) {
      shinyjs::disable("downloadGeneList")
    } else {
      shinyjs::enable("downloadGeneList")
    }
  })

  #create custom name for the results
  customName <- reactive(paste(input$selectDiffexCondition, input$selectDiffex, sep = "_"))
  # Download the differential expression results table
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste(customName(), Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      results <- vals$diffexgenelist
      colnames(results) <- paste(customName(), colnames(results), sep = "_")
      utils::write.csv(results, file)
    }
  )

  #save results wrt to custom name
  observeEvent(input$saveResults, {
    if (input$ResultsName == ""){
      shinyalert::shinyalert("Error!", "Specify name of the results.", type = "error")
    } else {
      withBusyIndicatorServer("saveResults", {
        ResultsName <- gsub(" ", "_", input$ResultsName)
        if (!is.null(diffExValues$diffExList)) {
          if (ResultsName %in% diffExValues$diffExList) {
            shinyalert::shinyalert("Error!", "name already exists. Please use a unique result name", type = "error")
          } else {
            diffExValues$index <- diffExValues$index + 1
            diffExValues$diffExList[diffExValues$index] <- ResultsName
            vals$counts <- saveDiffExResults(inSCE = vals$counts,
                                             diffex = vals$diffexgenelist,
                                             name = input$ResultsName,
                                             method = input$selectDiffex)
          }
        } else {
          diffExValues$index <- diffExValues$index + 1
          diffExValues$diffExList[diffExValues$index] <- ResultsName
          vals$counts <- saveDiffExResults(inSCE = vals$counts,
                                           diffex = vals$diffexgenelist,
                                           name = input$ResultsName,
                                           method = input$selectDiffex)
        }
      })
    }
  })

  #dynamically create a list of names of the results
  output$savedRes <- renderUI({
    if (!is.null(vals$counts)) {
      if (is.null(diffExValues$diffExList)) {
        savedObjResults <- gsub("_padj$", "", colnames(rowData(vals$counts))[grepl("_padj$", colnames(rowData(vals$counts)))])
        diffExValues$index <- length(savedObjResults)
        diffExValues$diffExList[seq_len(diffExValues$index)] <- savedObjResults[seq_len(diffExValues$index)]
      }
      selectizeInput("savedDiffExResults", "Select available results",
                     choices = diffExValues$diffExList)
    }
  })

  output$saveDiffResultsNote <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Note: Use a unique name to save results each time.")))
  })

  #load specific result according to users' input
  observeEvent(input$loadResults, {
    if (!is.null(input$savedDiffExResults)) {
      df <- data.frame(rowData(vals$counts))
      #extract all columns except biomarker columns by checking for unique values == 2
      listColNames <- names(which(apply(df, 2, function(a) length(unique(a)) == 2) == FALSE))
      #arrange your df according to these extracted names
      df <- df[, listColNames]
      #find all columns matching with the user's input.
      #sub() here is used to extract columns with user defined inputs
      df <- df[, which(sub("_[^_]+$", "", listColNames) == input$savedDiffExResults)]
      #remove the saved results names from columns
      colnames(df) <- gsub(paste0(input$savedDiffExResults, "_"), "", colnames(df))
      filterCol <- colnames(df)[grepl("padj$", colnames(df))]
      #order the padj column to get the top significant genes
      orderedRows <- rownames(df)[order(df[, filterCol])[seq_len(nrow(df))]]
      df <- df[orderedRows, ]
      vals$diffexgenelist <- df
      vals$diffexheatmapplot <- NULL
    }
  })

  output$BioNgenes <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Max genes: "), nrow(vals$diffexgenelist), sep = ""))
  })

  output$logFCBioRange <- renderUI({
    req(vals$diffexgenelist)
    if (input$selectDiffex != 'ANOVA') {
      logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
      if (length(logFCIndex)) {
        #for DESeq2 with more than 1 covariate, choose the first column
        logFCIndex <- 1
      }
      minlogFC <- paste(em("Min logFC : "), round(min(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      maxlogFC <- paste(em("Max logFC : "), round(max(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      HTML(paste(minlogFC, maxlogFC, sep = '<br/>'))
    }
  })

  #save biomarker in rowData() wrt name and conditions.
  observeEvent(input$saveBiomarker, {
    if (input$biomarkerName == ""){
      shinyalert::shinyalert("Error!", "Specify biomarker name.", type = "error")
    } else {
      withBusyIndicatorServer("saveBiomarker", {
        req(vals$diffexgenelist)
        biomarkerName <- gsub(" ", "_", input$biomarkerName)
        if (anyDuplicated(biomarkerName)) {
          shinyalert::shinyalert("Error", "name already exists. Please use a unique result name",
                                 type = "error")
        }
        if (input$applyAbslogFC == TRUE) {
          absLogFC <- abs(input$selectlogFC)
        } else {
          absLogFC <- input$selectlogFC
        }
        if (input$selectBioNGenes > nrow(vals$diffexgenelist)) {
          stop("Max value exceeded for Input.")
        }
        if (input$applyBioCutoff1 == TRUE & input$applyBioCutoff2 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = absLogFC,
                                          pVal = input$selectAdjPVal)
        } else if (input$applyBioCutoff1 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = NULL,
                                          pVal = input$selectAdjPVal)
        } else if (input$applyBioCutoff2 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = absLogFC,
                                          pVal = NULL)
        } else {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = input$biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = NULL,
                                          pVal = NULL)
        }
        vals$diffexBmName <- TRUE
      })
    }
  })

  observe({
    output$bioMarkerNote <- renderUI({
      req(vals$counts)
      req(isolate(input$biomarkerName))
      if (vals$diffexBmName) {
        isolate({
          biomarkerName <- gsub(" ", "_", input$biomarkerName)
          countBioGenes <- count(rowData(vals$counts)[, biomarkerName] == 1)
          HTML(paste("Saved ", countBioGenes, " genes after applying the selected filter(s)", sep = ""))
        })
      }
    })
  })

  #-----------------------------------------------------------------------------
  # Page 5.2: MAST ####
  #-----------------------------------------------------------------------------
  ## MAST - condition determination method1 ####
  output$mastC1G1UI <- renderUI({
    if(!is.null(vals$counts) &
       !input$mastC1Class == "None"){
      classCol <- colData(vals$counts)[[input$mastC1Class]]
      classChoices <- sort(as.vector(unique(classCol)))
      selectInput(inputId = "mastC1G1", label = "Select Condition(s)",
                  choices = classChoices, multiple = TRUE)
    } else {
      selectInput(inputId = "mastC1G1", label = "Select Condition(s)",
        choices = NULL, multiple = TRUE)
    }
  })

  output$mastC1G2UI <- renderUI({
    if(!is.null(vals$counts) &
        !input$mastC1Class == "None"){
      classCol <- colData(vals$counts)[[input$mastC1Class]]
      classChoices <- sort(as.vector(unique(classCol)))
      selectInput(inputId = "mastC1G2", label = "Select Condition(s)",
        choices = classChoices, multiple = TRUE)
    } else {
      selectInput(inputId = "mastC1G2", label = "Select Condition(s)",
        choices = NULL, multiple = TRUE)
    }
  })

  output$mastC1G1CellCheckUI <- renderUI({
    if(!is.null(input$mastC1G1) &
       length(input$mastC1G1) > 0){
      g1Idx <- colData(vals$counts)[[input$mastC1Class]] %in% input$mastC1G1
      g1Cells <- colnames(vals$counts)[g1Idx]
      g1CellsText <- paste(g1Cells, collapse = "\n")
      textAreaInput("mastC1G1CellCheck", "Cells selected:", g1CellsText,
        height = '100px', placeholder = "Nothing selected")
    } else {
      textAreaInput("mastC1G1CellCheck", "Cells selected:", NULL,
        height = '100px', placeholder = "Nothing selected")
    }
  })

  output$mastC1G2CellCheckUI <- renderUI({
    if(!is.null(input$mastC1G2) &
        length(input$mastC1G2) > 0){
      g2Idx <- colData(vals$counts)[[input$mastC1Class]] %in% input$mastC1G2
      g2Cells <- colnames(vals$counts)[g2Idx]
      g2CellsText <- paste(g2Cells, collapse = "\n")
      textAreaInput("mastC1G2CellCheck", "Cells selected:", g2CellsText,
        height = '100px', placeholder = "Nothing selected")
    } else {
      textAreaInput("mastC1G2CellCheck", "Cells selected:", NULL,
        height = '100px', placeholder = "Nothing selected")
    }
  })

  output$mastC1G1NCell <- renderUI({
    if(!is.null(input$mastC1G1CellCheck)){
      if(!input$mastC1G1CellCheck == ""){
        cellList <- str_trim(scan(text = input$mastC1G1CellCheck,
          sep='\n', what = 'character'))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " cell(s) selected.")
    span(msg, style = 'margin-left:10px')
  })

  output$mastC1G2NCell <- renderUI({
    if(!is.null(input$mastC1G2CellCheck)){
      if(!input$mastC1G2CellCheck == ""){
        cellList <- str_trim(scan(text = input$mastC1G2CellCheck,
          sep='\n', what = 'character'))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " cell(s) selected.")
    span(msg, style = 'margin-left:10px')
  })
  ## MAST - condition determination method2 ####
  ## condition 1 table operation vvvv
  output$mastC2G1Table <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      data.frame(colData(vals$counts)[,input$mastC2G1Col])
    }
  }, filter = "top", server = TRUE
  )
  mastC2G1Table_proxy <- DT::dataTableProxy("mastC2G1Table")

  observeEvent(input$mastC2G1Table_addAll, {
    DT::selectRows(mastC2G1Table_proxy,
      sort(unique(c(input$mastC2G1Table_rows_selected,
                    input$mastC2G1Table_rows_all))))
  })

  observeEvent(input$mastC2G1Table_clear, {
    DT::selectRows(mastC2G1Table_proxy, NULL)
  })

  output$mastC2G1info <- renderUI({
    nCell <- length(input$mastC2G1Table_rows_selected)
    p(paste0("Totally ", nCell, " cells selected for ", input$mastG1Name))
  })
  ## condition 1 table operation ^^^^
  ## condition 2 table operation vvvv
  output$mastC2G2Table <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      data.frame(colData(vals$counts)[,input$mastC2G2Col])
    }
  }, filter = "top", server = TRUE
  )
  mastC2G2Table_proxy <- DT::dataTableProxy("mastC2G2Table")

  observeEvent(input$mastC2G2Table_addAll, {
    DT::selectRows(mastC2G2Table_proxy,
      sort(unique(c(input$mastC2G2Table_rows_selected,
        input$mastC2G2Table_rows_all))))
  })

  observeEvent(input$mastC2G2Table_clear, {
    DT::selectRows(mastC2G2Table_proxy, NULL)
  })

  output$mastC2G2info <- renderUI({
    nCell <- length(input$mastC2G2Table_rows_selected)
    p(paste0("Totally ", nCell, " cells selected for ", input$mastG2Name))
  })
  ## condition 2 table operation ^^^^
  ## MAST - condition determination method3 ####
  output$mastC3G1NCell <- renderUI({
    if(!is.null(input$mastC3G1Cell)){
      if(!input$mastC3G1Cell == ""){
        cellList <- str_trim(scan(text = input$mastC3G1Cell,
          sep='\n', what = 'character'))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " valid cell name(s) entered.")
    span(msg, style = 'margin-left:10px')
  })

  output$mastC3G2NCell <- renderUI({
    if(!is.null(input$mastC3G2Cell)){
      if(!input$mastC3G2Cell == ""){
        cellList <- str_trim(scan(text = input$mastC3G2Cell,
          sep='\n', what = 'character'))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " valid cell name(s) entered.")
    span(msg, style = 'margin-left:10px')
  })
  ## MAST - other check ####
  output$mastCompNameUI <- renderUI({
    if(!is.null(vals$counts)){
      if("MAST" %in% names(S4Vectors::metadata(vals$counts))){
        nRes <- length(names(S4Vectors::metadata(vals$counts)$MAST))
        autoinc <- nRes + 1
      } else {
        autoinc <- 1
      }
    } else {
      autoinc <- 1
    }
    textInput("mastCompName", "Experiment Name:", paste0("Comparison", autoinc),
              placeholder = 'Required.')
  })

  mastNameIsDup <- reactive({
    if(!is.null(vals$counts) &&
        "MAST" %in% names(metadata(vals$counts)) &&
        !is.null(input$mastCompName)){
      allRes <- names(metadata(vals$counts)$MAST)
      input$mastCompName %in% allRes
    } else {
      FALSE
    }
  })

  runMASTfromShiny <- function(){
    withBusyIndicatorServer("runMAST", {
      if(input$mastCondMethod == 1){
        vals$counts <- runMAST(inSCE = vals$counts,
                               useAssay = input$mastAssay, class = input$mastC1Class,
                               classGroup1 = input$mastC1G1, classGroup2 = input$mastC1G2,
                               groupName1 = input$mastG1Name, groupName2 = input$mastG2Name,
                               comparisonName = input$mastCompName,
                               useThresh = input$useAdaptThresh, freqExpressed = input$mastFreq,
                               log2fcThreshold = input$mastFCThresh,
                               fdrThreshold = input$mastFDRThresh, onlyPos = input$mastPosOnly)
      } else if(input$mastCondMethod == 2){
        vals$counts <- runMAST(inSCE = vals$counts,
                               useAssay = input$mastAssay,
                               index1 = input$mastC2G1Table_rows_selected,
                               index2 = input$mastC2G2Table_rows_selected,
                               groupName1 = input$mastG1Name, groupName2 = input$mastG2Name,
                               comparisonName = input$mastCompName,
                               useThresh = input$useAdaptThresh, freqExpressed = input$mastFreq,
                               log2fcThreshold = input$mastFCThresh,
                               fdrThreshold = input$mastFDRThresh, onlyPos = input$mastPosOnly)
      } else {
        g1CellList <- str_trim(scan(text = input$mastC3G1Cell,
                                    sep='\n', what = 'character'))
        g1CellList <- sort(unique(g1CellList))
        g2CellList <- str_trim(scan(text = input$mastC3G2Cell,
                                    sep='\n', what = 'character'))
        g2CellList <- sort(unique(g2CellList))
        vals$counts <- runMAST(inSCE = vals$counts,
                               useAssay = input$mastAssay, index1 = g1CellList,
                               index2 = g2CellList, groupName1 = input$mastG1Name,
                               groupName2 = input$mastG2Name, comparisonName = input$mastCompName,
                               useThresh = input$useAdaptThresh, freqExpressed = input$mastFreq,
                               log2fcThreshold = input$mastFCThresh,
                               fdrThreshold = input$mastFDRThresh, onlyPos = input$mastPosOnly)
      }
      shinyalert::shinyalert("Success",
                             text = "MAST Differential Expression completed.", type = "success")
      allResName <- names(metadata(vals$counts)$MAST)
      updateSelectInput(session, "mastResSel", choices = allResName)
    })
  }
  ## MAST - apply calculation ####
  observeEvent(input$runMAST, {
    if (is.null(vals$counts)){
      shinyalert("Error!", "Upload data first.", type = "error")
    } else if(input$mastCompName == ""){
      shinyalert("Error!", "Please enter differential expression analysis name.", type = "error")
    } else {
      allRes <- names(metadata(vals$counts)$MAST)
      if(input$mastCompName %in% allRes){
        shinyalert("Warning",
                   "Entered differential experiment analysis name is already there.",
                   "warning", showCancelButton = TRUE,
                   confirmButtonText = "Overwrite",
                   callbackR = function(x){if(isTRUE(x)){runMASTfromShiny()}})
      } else {
          runMASTfromShiny()
      }
    }
  })

  output$mastResSelUI <- renderUI({
    if(!is.null(vals$counts)){
      res <- names(metadata(vals$counts)$MAST)
      selectInput("mastResSel", "Select Differential Expression Analysis", res)
    }
  })

  output$threshplot <- renderPlot({
    if(!is.null(vals$counts)){
      vals$thres <- thresholdGenes(inSCE = vals$counts,
        useAssay = input$mastAssay)
      par(mfrow = c(5, 4))
      plot(vals$thres)
      par(mfrow = c(1, 1))
    }
  }, height = 800)

  output$mastresults <- DT::renderDataTable({
    if(!is.null(input$mastResSel)){
      metadata(vals$counts)$MAST[[input$mastResSel]]$result
    }
  }, filter = 'top')

  isMastResult <- reactive(is.null(input$mastResSel) ||
      input$mastResSel == "")
  observe({
    if (isMastResult()) {
      shinyjs::disable("mastDownload")
    } else {
      shinyjs::enable("mastDownload")
    }
  })

  output$mastDownload <- downloadHandler(
    filename = function() {
      paste0("mastResult_", input$mastResSel, ".csv")
    },
    content = function(file) {
      fullTable <- metadata(vals$counts)$MAST[[input$mastResSel]]$result
      filteredTable <- fullTable[input$mastresults_rows_all,]
      utils::write.csv(filteredTable, file, row.names = FALSE)
    }
  )

  output$mastVioTotalUI <- renderUI({
    topN <- input$mastVioNrow * input$mastVioNcol
    p(as.character(topN))
  })

  output$hurdleviolin <- renderPlot({
    if(!is.null(input$mastResSel) &&
       !input$mastResSel == ""){
      plotMASTViolin(inSCE = vals$counts, useResult = input$mastResSel,
                     threshP = input$mastVioUseThresh,
                     nrow = input$mastVioNrow, ncol = input$mastVioNcol)
    }
  })

  output$mastRegTotalUI <- renderUI({
    topN <- input$mastRegNrow * input$mastRegNcol
    p(as.character(topN))
  })

  output$hurdlelm <- renderPlot({
    if(!is.null(input$mastResSel) &&
       !input$mastResSel == ""){
      plotMASTRegression(inSCE = vals$counts, useResult = input$mastResSel,
                         threshP = input$mastRegUseThresh,
                         nrow = input$mastRegNrow, ncol = input$mastRegNcol)
    }
  })

  output$hurdleHeatmap <- renderPlot({
    if(!is.null(input$mastResSel) &&
       !input$mastResSel == ""){
      if(length(input$mastHMSplitCol) == 1 && input$mastHMSplitCol == 'None'){
        colSplitBy <- NULL
      } else if(!'None' %in% input$mastHMSplitCol){
        colSplitBy <- input$mastHMSplitCol
      } else {
        whereNone <- input$mastHMSplitCol %in% 'None'
        colSplitBy <- input$mastHMSplitCol
        colSplitBy <- colSplitBy[!whereNone]
      }
      if(length(input$mastHMSplitRow) == 1 && input$mastHMSplitRow == 'None'){
        rowSplitBy <- NULL
      } else if(!'None' %in% input$mastHMSplitRow){
        rowSplitBy <- input$mastHMSplitRow
      } else {
        whereNone <- input$mastHMSplitRow %in% 'None'
        rowSplitBy <- input$mastHMSplitRow
        rowSplitBy <- rowSplitBy[!whereNone]
      }
      plotMASTHeatmap(inSCE = vals$counts, useResult = input$mastResSel,
        onlyPos = input$mastHMPosOnly, log2fcThreshold = input$mastHMFC,
        fdrThreshold = input$mastHMFDR, rowDataName = input$mastHMrowData,
        colDataName = input$mastHMcolData, colSplitBy = colSplitBy,
        rowSplitBy = rowSplitBy)
    }
  })

  #-----------------------------------------------------------------------------
  # Page 5.3: MAST - Find Marker ####
  #-----------------------------------------------------------------------------
  # MAST findMarker RUN ####
  observeEvent(input$runMASTFM, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("runMASTFM", {
        vals$counts <- findMarkerDiffExp(vals$counts,
          useAssay = input$mastFMAssay, cluster = input$mastFMCluster,
          log2fcThreshold = input$mastFMLogFC, fdrThreshold = input$mastFMFDR,
          useThresh = input$mastFMUseThresh, freqExpressed = input$mastFMFreq)
        shinyalert::shinyalert("Success", "MAST Find Marker completed.",
                               "success")
      })
    }
  })
  # MAST findMarker ResultTable ####
  output$mastFMResClusterUI <- renderUI({
    if(!is.null(vals$counts) &&
       !is.null(input$mastFMCluster) &&
       'findMarker' %in% names(metadata(vals$counts))){
      allCluster <- colData(vals$counts)[[input$mastFMCluster]]
      allCluster <- sort(as.vector(unique(allCluster)))
      checkboxGroupInput("mastFMResCluster", "View Results for Cluster",
                         choices = c('All', allCluster), selected = 'All',
                         inline = TRUE)
    }
  })

  output$mastFMResTable <- DT::renderDataTable({
    if(!is.null(vals$counts) &&
       'findMarker' %in% names(metadata(vals$counts))){
      fullTable <- metadata(vals$counts)$findMarker
      if('All' %in% input$mastFMResCluster){
        fullTable
      } else {
        fullTable[fullTable[[input$mastFMCluster]] %in% input$mastFMResCluster,]
      }
    }
  }, filter = "top")

  isMastFMResult <- reactive(!is.null(vals$counts) &&
                             !is.null(metadata(vals$counts)$findMarker))
  observe({
    if (isMastFMResult()) {
      shinyjs::enable("mastFMDownload")
    } else {
      shinyjs::disable("mastFMDownload")
    }
  })

  output$mastFMDownload <- downloadHandler(
    filename = function() {
      paste0("mastMarkerResult_", input$mastFMCluster, ".csv")
    },
    content = function(file) {
      fullTable <- metadata(vals$counts)$findMarker
      filteredTable <- fullTable[input$mastFMResTable_rows_all,]
      utils::write.csv(filteredTable, file, row.names = FALSE)
    }
  )

  # MAST findMarker Heatmap ####
  output$mastFMHMAssayUI <- renderUI({
    if(!is.null(vals$counts)){
      allAssay <- assayNames(vals$counts)
      selectInput('mastFMHMAssay', "Assay to plot", allAssay,
                  selected = input$mastFMAssay)
    }
  })

  output$mastFMHeatmap <- renderPlot({
    if(!is.null(vals$counts) &&
       'findMarker' %in% names(metadata(vals$counts)) &&
       !is.null(input$mastFMHMAssay)){
      plotMarkerDiffExp(inSCE = vals$counts, useAssay = input$mastFMHMAssay,
        orderBy = input$mastFMHMOrder, log2fcThreshold = input$mastFMHMFC,
        fdrThreshold = input$mastFMHMFDR, decreasing = input$mastFMHMdec,
        rowDataName = input$mastFMHMrowData, colDataName = input$mastFMHMcolData
      )
    }
  })
  #-----------------------------------------------------------------------------
  # Page 6: Pathway Activity Analysis
  #-----------------------------------------------------------------------------

  output$selectPathwayGeneLists <- renderUI({
    if (input$genelistSource == "Manual Input"){
      if (!is.null(vals$counts)){
        #fn to check if each column is 1 and 0 only
        biomarkercols <- names(which(apply(rowData(vals$counts), 2, function(a) length(unique(a)) == 2) == TRUE))
        selectizeInput("pathwayGeneLists", "Select Gene List(s):",
                       biomarkercols, multiple = TRUE)
      } else {
        h4("Note: upload data.")
      }
    } else {
      selectInput("pathwayGeneLists", "Select Gene List(s):",
                  c("ALL", names(c2BroadSets)), multiple = TRUE)
    }
  })

  output$selectNumTopPaths <- renderUI({
    if (!is.null(input$pathwayGeneLists)) {
      if ("ALL" %in% input$pathwayGeneLists & input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
        sliderInput("pickNtopPaths", "Number of top pathways:", min = 5,
                    max = length(c2BroadSets), value = 25, step = 5)
      }
    }
  })

  observeEvent(input$pathwayRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("pathwayRun", {
        vals$gsvaRes <- gsvaSCE(inSCE = vals$counts,
                                useAssay = input$pathwayAssay,
                                pathwaySource = input$genelistSource,
                                pathwayNames = input$pathwayGeneLists)
      })
    }
  })

  observe({
    if (length(input$pathwayPlotVar) == 1 & !(is.null(vals$gsvaRes))){
      fit <- limma::lmFit(vals$gsvaRes, stats::model.matrix(~factor(colData(vals$counts)[, input$pathwayPlotVar])))
      fit <- limma::eBayes(fit)
      toptableres <- limma::topTable(fit, number = nrow(vals$gsvaRes))
      temptable <- cbind(rownames(toptableres), toptableres)
      rownames(temptable) <- NULL
      colnames(temptable)[1] <- "Pathway"
      vals$gsvaLimma <- temptable
    } else {
      vals$gsvaLimma <- NULL
    }
  })

  output$pathwaytable <- DT::renderDataTable({
    if (!is.null(vals$gsvaLimma)){
      if (!is.null(input$pathwayGeneLists) & "ALL" %in% input$pathwayGeneLists & input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
        vals$gsvaLimma[1:min(input$pickNtopPaths, nrow(vals$gsvaLimma)), , drop = FALSE]
      } else {
        vals$gsvaLimma
      }
    }
  }, options = list(scrollX = TRUE, pageLength = 30))

  output$pathwayPlot <- renderPlot({
    if (!(is.null(vals$gsvaRes))){
      if (input$genelistSource == "MSigDB c2 (Human, Entrez ID only)" & "ALL" %in% input$pathwayGeneLists & !(is.null(vals$gsvaLimma))){
        tempgsvares <- vals$gsvaRes[as.character(vals$gsvaLimma$Pathway[1:min(input$pickNtopPaths, nrow(vals$gsvaLimma))]), , drop = FALSE]
      } else if (input$genelistSource == "MSigDB c2 (Human, Entrez ID only)" & !("ALL" %in% input$pathwayGeneLists)) {
        tempgsvares <- vals$gsvaRes
      } else {
        tempgsvares <- vals$gsvaRes[1:input$pickNtopPaths, , drop = FALSE]
      }
      if (input$pathwayOutPlot == "Violin" & length(input$pathwayPlotVar) > 0){
        tempgsvares <- tempgsvares[1:min(49, input$pickNtopPaths, nrow(tempgsvares)), , drop = FALSE]
        gsvaPlot(inSCE = vals$counts,
                 gsvaData = tempgsvares,
                 plotType = input$pathwayOutPlot,
                 condition = input$pathwayPlotVar)
      } else if (input$pathwayOutPlot == "Heatmap"){
        gsvaPlot(inSCE = vals$counts,
                 gsvaData = tempgsvares,
                 plotType = input$pathwayOutPlot,
                 condition = input$pathwayPlotVar)
      }
    }
  })

  #save pathawy activity results in the colData
  observeEvent(input$savePathway, {
    if (!(is.null(vals$gsvaRes))){
      if (all(colnames(vals$counts) == colnames(vals$gsvaRes))){
        #if we have limma results
        if (!(is.null(vals$gsvaLimma))){
          tempdf <- DataFrame(t(vals$gsvaRes[vals$gsvaLimma$Pathway[1:input$pickNtopPaths], , drop = FALSE]))
        } else {
          tempdf <- DataFrame(t(vals$gsvaRes[1:input$pickNtopPaths, , drop = FALSE]))
        }
        tempdf <- tempdf[, !(colnames(tempdf) %in% colnames(colData(vals$counts))), drop = FALSE]
        colData(vals$counts) <- cbind(colData(vals$counts), tempdf)
        updateColDataNames()
      }
    } else {
      shinyalert::shinyalert("Error!", "Run pathway first.", type = "error")
    }
  })

  #disable downloadPathway button if the pathway data doesn't exist
  isPathwayResult <- reactive(is.null(vals$gsvaRes))
  observe({
    if (isPathwayResult()) {
      shinyjs::disable("downloadPathway")
    } else {
      shinyjs::enable("downloadPathway")
    }
  })

  #download mast results
  output$downloadPathway <- downloadHandler(
    filename = function() {
      paste("pathway_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(vals$gsvaRes, file)
    }
  )

  #-----------------------------------------------------------------------------
  # Page 6.2 : Enrichment Analysis - EnrichR
  #-----------------------------------------------------------------------------

  enrichRfile <- reactive(read.csv(input$enrFile$datapath,
                                   header = input$header,
                                   sep = input$sep,
                                   quote = input$quote,
                                   row.names = 1))
  output$enrBioGenes <- renderUI({
    if (!is.null(vals$counts)) {
      selectInput("selEnrBioGenes", "Select Gene List(s):",
                  names(which(apply(rowData(vals$counts), 2, function(a) length(unique(a)) == 2) == TRUE)))
    }
  })
  dbs <- reactive({
    if (internetConnection){
      enrDatabases <- enrichR::listEnrichrDbs()$libraryName
    } else {
      enrDatabases <- ""
    }
    if (is.null(input$enrichDb)){
      dbs <- enrDatabases
    } else {
      if (any(input$enrichDb %in% "ALL")){
        dbs <- enrDatabases
      } else {
        dbs <- input$enrichDb
      }
    }
  })

  #count_db <- reactive(length(dbs()))
  observeEvent (input$enrichRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer ("enrichRun", {
        tryCatch ({
          if (input$geneListChoice == "selectGenes"){
            genes <- input$enrichGenes
          } else if (input$geneListChoice == "geneFile"){
            req(input$enrFile)
            genes <- rownames(enrichRfile())
          } else  {
            genes <- rownames(vals$counts)[SingleCellExperiment::rowData(vals$counts)[, input$selEnrBioGenes] == 1]
          }
          vals$enrichRes <- enrichRSCE(inSCE = vals$counts,
                                       glist = genes,
                                       db = dbs())
        }, error = function(e){
          shinyalert::shinyalert("Error!", e$message, type = "error")
        })
      })
    }
  })

  output$enrTabs <- renderUI({
    req(vals$enrichRes)
    isoDbs <- isolate(dbs())
    nTabs <- length(isoDbs)
    #create tabPanel with datatable in it
    myTabs <- lapply(seq_len((nTabs)), function(i) {
      tabPanel(paste0(isoDbs[i]),
               DT::dataTableOutput(paste0(isoDbs[i]))
      )
    })
    do.call(tabsetPanel, myTabs)
  })

  #create datatables
  observe({
    req(vals$enrichRes)
    isoDbs <- isolate(dbs())
    enrResults <- vals$enrichRes[, c(1:10)] %>%
      mutate(Database_selected =
               paste0("<a href='", vals$enrichRes[, 11],
                      "' target='_blank'>",
                      vals$enrichRes[, 1], "</a>"))
    lapply(seq_len(length(isoDbs)), function(i){
      output[[paste0(isoDbs[i])]] <- DT::renderDataTable({
        DT::datatable({
          enr <- enrResults[which(vals$enrichRes[, 1] %in% isoDbs[i]), ]
        }, escape = FALSE, options = list(scrollX = TRUE, pageLength = 30), rownames = FALSE)
      })
    })
  })

  #disable the downloadEnrichR button if the result doesn't exist
  isResult <- reactive(is.null(vals$enrichRes))
  observe({
    if (isResult()) {
      shinyjs::disable("downloadEnrichR")
    } else {
      shinyjs::enable("downloadEnrichR")
    }
  })

  output$downloadEnrichR <- downloadHandler(
    filename = function() {
      paste("enrichR-results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(vals$enrichRes, file)
    },
    contentType = "text/csv"
  )

  #-----------------------------------------------------------------------------
  # Page 7: Subsampling
  #-----------------------------------------------------------------------------

  #Run subsampling analysis
  observeEvent(input$runSubsampleDepth, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSubsampleDepth", {
        vals$subDepth <- DownsampleDepth(originalData = vals$counts,
                                         useAssay = input$depthAssay,
                                         minCount = input$minCount,
                                         minCells = input$minCells,
                                         maxDepth = 10 ^ input$maxDepth,
                                         realLabels = input$selectReadDepthCondition,
                                         depthResolution = input$depthResolution,
                                         iterations = input$iterations)

        output$depthDone <- renderPlot({
          plot(apply(vals$subDepth[, , 1], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Number of detected genes",
               main = "Number of dected genes by sequencing depth")
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$minEffectDone <- renderPlot({
          plot(apply(vals$subDepth[, , 2], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Average significant effect size",
               ylim = c(0, 2))
          lines(apply(vals$subDepth[, , 2], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 2], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$sigNumDone <- renderPlot({
          plot(apply(vals$subDepth[, , 3], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Number of significantly DiffEx genes")
          lines(apply(vals$subDepth[, , 3], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 3], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
      })
    }
  })

  observeEvent(input$runSubsampleCells, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSubsampleCells", {
        if (input$useReadCount){
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           useAssay = input$cellsAssay,
                                           realLabels = input$selectCellNumCondition,
                                           totalReads = sum(SummarizedExperiment::assay(vals$counts, input$cellsAssay)),
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        else{
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           useAssay = input$cellsAssay,
                                           realLabels = input$selectCellNumCondition,
                                           totalReads = input$totalReads,
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        output$cellsDone <- renderPlot({
          plot(apply(vals$subCells[, , 1], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of virtual cells", ylab = "Number of detected genes",
               main = "Number of dected genes by cell number")
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$minEffectCells <- renderPlot({
          plot(apply(vals$subCells[, , 2], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of virtual cells", ylab = "Average significant effect size",
               ylim = c(0, 2))
          lines(apply(vals$subCells[, , 2], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 2], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$sigNumCells <- renderPlot({
          plot(apply(vals$subCells[, , 3], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of vitual cells", ylab = "Number of significantly DiffEx genes")
          lines(apply(vals$subCells[, , 3], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 3], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
      })
    }
  })

  #Run differential power analysis
  observeEvent(input$runSnapshot, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSnapshot", {
        vals$snapshot <- iterateSimulations(originalData = vals$counts,
                                            useAssay = input$snapshotAssay,
                                            realLabels = input$selectSnapshotCondition,
                                            totalReads = input$numReadsSnap,
                                            cells = input$numCellsSnap,
                                            iterations = input$iterationsSnap)
        vals$effectSizes <- calcEffectSizes(countMatrix = SummarizedExperiment::assay(vals$counts, input$snapshotAssay), condition = colData(vals$counts)[, input$selectSnapshotCondition])
        output$Snaplot <- renderPlot({
          plot(apply(vals$snapshot, 1, function(x){sum(x <= 0.05) / length(x)}) ~ vals$effectSizes,
               xlab = "Cohen's d effect size", ylab = "Detection power", lwd = 4, main = "Power to detect diffex by effect size")
        })
      })
    }
  })

    #-----------------------------------------------------------------------------
    # Page 8: Seurat Workflow
    #-----------------------------------------------------------------------------

  #Perform normalization
  observeEvent(input$normalize_button, {
    req(vals$counts)
    withProgress(message = "Normalizing", max = 1, value = 1, {
      vals$counts <- seuratNormalizeData(inSCE = vals$counts,
                                         useAssay = input$seuratSelectNormalizationAssay,
                                         normAssayName = "seuratNormData",
                                         normalizationMethod = input$normalization_method,
                                         scaleFactor = as.numeric(input$scale_factor))
      updateAssayInputs()
      vals$counts <- .seuratInvalidate(inSCE = vals$counts)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Normalize Data" = "danger"))
    shinyjs::enable(selector = "div[value='Scale Data']")
    showNotification("Normalization Complete")
  })

  #Perform scaling
  observeEvent(input$scale_button, {
    req(vals$counts)
    withProgress(message = "Scaling", max = 1, value = 1, {
      vals$counts <- seuratScaleData(inSCE = vals$counts,
                                     useAssay = "seuratNormData",
                                     scaledAssayName = "seuratScaledData",
                                     model = input$model.use,
                                     scale = input$do.scale,
                                     center = input$do.center,
                                     scaleMax = input$scale.max)
      updateAssayInputs()
      vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "danger"))
    shinyjs::enable(selector = "div[value='Highly Variable Genes']")
    showNotification("Scale Complete")
  })

  #Find HVG
  observeEvent(input$find_hvg_button, {
    req(vals$counts)
    withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
      vals$counts <- seuratFindHVG(inSCE = vals$counts,
                                   useAssay = "seuratScaledData",
                                   hvgMethod = input$hvg_method,
                                   hvgNumber = as.numeric(input$hvg_no_features))

      vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE)
    })
    withProgress(message = "Plotting HVG", max = 1, value = 1, {
      output$plot_hvg <- renderPlotly({
        plotly::ggplotly(seuratPlotHVG(vals$counts))
      })
    })
    updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "danger"))
    shinyjs::enable(selector = "div[value='Dimensionality Reduction']")
    showNotification("Find HVG Complete")
  })

  #Display highly variable genes
  output$hvg_output <- renderText({
    if (!is.null(vals$counts)) {
      if (!is.null(vals$counts@metadata$seurat$obj)) {
        if (length(slot(vals$counts@metadata$seurat$obj, "assays")[["RNA"]]@var.features) > 0) {
          .seuratGetVariableFeatures(vals$counts, input$hvg_no_features_view)
        }
      }
    }
  })

  #Run PCA
  observeEvent(input$run_pca_button, {
    req(vals$counts)
    withProgress(message = "Running PCA", max = 1, value = 1, {
      vals$counts <- seuratPCA(inSCE = vals$counts,
                               useAssay = "seuratScaledData",
                               reducedDimName = "seuratPCA",
                               nPCs = input$pca_no_components)

      vals$counts@metadata$seurat$count_pc <- dim(convertSCEToSeurat(vals$counts)[["pca"]])[2]
      vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE)
    })
    withProgress(message = "Plotting PCA", max = 1, value = 1, {
      output$plot_pca <- renderPlotly({
        plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                            useReduction = "pca",
                            showLegend = FALSE))
      })
    })
    if (input$pca_compute_elbow) {
      withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
        updateSliderInput(session = session, inputId = "pca_significant_pc_slider", value = .computeSignificantPC(vals$counts))
        output$plot_elbow_pca <- renderPlotly({
          seuratElbowPlot(inSCE = vals$counts,
                          significantPC = .computeSignificantPC(vals$counts))
        })
        output$pca_significant_pc_output <- renderText({
          .computeSignificantPC(vals$counts)
        })
      })
    }
    if (input$pca_compute_jackstraw) {
      withProgress(message = "Generating JackStraw Plot", max = 1, value = 1, {
        vals$counts <- seuratComputeJackStraw(inSCE = vals$counts,
                                              useAssay = "seuratScaledData",
                                              dims = input$pca_no_components)
        output$plot_jackstraw_pca <- renderPlotly({
          plotly::ggplotly(seuratJackStrawPlot(inSCE = vals$counts,
                              dims = input$pca_no_components))
        })
      })
    }
    if (input$pca_compute_heatmap) {
      withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
        vals$counts@metadata$seurat$heatmap_pca <- seuratComputeHeatmap(inSCE = vals$counts,
                                                                        useAssay = "seuratScaledData",
                                                                        useReduction = "pca",
                                                                        dims = input$pca_no_components,
                                                                        nfeatures = input$pca_compute_heatmap_nfeatures,
                                                                        combine = FALSE,
                                                                        fast = FALSE)
        output$plot_heatmap_pca <- renderPlot({
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                            dims = input$pca_no_components,
                            ncol = 2,
                            labels = c("PC1", "PC2", "PC3", "PC4"))
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_pca", choices = .getComponentNames(vals$counts@metadata$seurat$count_pc, "PC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))

    #Enable/Disable PCA plot panels not selected for computation (ElbowPlot, JackStraw or Heatmap)
    shinyjs::enable(
      selector = ".seurat_pca_plots a[data-value='PCA Plot']")

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='Elbow Plot']",
      condition = input$pca_compute_elbow)

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='JackStraw Plot']",
      condition = input$pca_compute_jackstraw)

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='Heatmap Plot']",
      condition = input$pca_compute_heatmap)

    shinyjs::enable(
      selector = "div[value='tSNE/UMAP']")

    showNotification("PCA Complete")
  })

  #Run ICA
  observeEvent(input$run_ica_button, {
    req(vals$counts)
    withProgress(message = "Running ICA", max = 1, value = 1, {
      vals$counts <- seuratICA(inSCE = vals$counts,
                               useAssay = "seuratScaledData",
                               nics = input$ica_no_components)

      vals$counts@metadata$seurat$count_ic <- dim(convertSCEToSeurat(vals$counts)[["ica"]])[2]
      vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE)
    })
    withProgress(message = "Plotting ICA", max = 1, value = 1, {
      output$plot_ica <- renderPlotly({
        plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                            useReduction = "ica",
                            showLegend = FALSE))
      })
    })
    if (input$ica_compute_heatmap) {
      withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
        vals$counts@metadata$seurat$heatmap_ica <- seuratComputeHeatmap(inSCE = vals$counts,
                                                                        useAssay = "seuratScaledData",
                                                                        useReduction = "ica",
                                                                        dims = input$ica_no_components,
                                                                        nfeatures = input$ica_compute_heatmap_nfeatures,
                                                                        combine = FALSE,
                                                                        fast = FALSE)
        output$plot_heatmap_ica <- renderPlot({
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                            dims = input$ica_no_components,
                            ncol = 2,
                            labels = c("IC1", "IC2", "IC3", "IC4"))
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_ica", choices = .getComponentNames(vals$counts@metadata$seurat$count_ic, "IC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))

    #Enable/Disable ICA plot panels not selected for computation (Heatmap)
    shinyjs::enable(
      selector = ".seurat_ica_plots a[data-value='ICA Plot']")

    shinyjs::toggleState(
      selector = ".seurat_ica_plots a[data-value='Heatmap Plot']",
      condition = input$ica_compute_heatmap)

    shinyjs::enable(
      selector = "div[value='tSNE/UMAP']")

    showNotification("ICA Complete")
  })

  #Find clusters
  observeEvent(input$find_clusters_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_clustering_method]])){
      withProgress(message = "Finding clusters", max = 1, value = 1, {
        vals$counts <- seuratFindClusters(inSCE = vals$counts,
                                          useAssay = "seuratScaledData",
                                          useReduction = input$reduction_clustering_method,
                                          dims = input$pca_significant_pc_slider,
                                          algorithm = input$algorithm.use,
                                          groupSingletons = input$group.singletons,
                                          resolution = input$resolution_clustering)
      })
      updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "danger"))
      showNotification("Find Clusters Complete")

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]])){
        withProgress(message = "Re-generating PCA plot with cluster labels", max = 1, value = 1,{
          output$plot_pca_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                useReduction = "pca",
                                showLegend = TRUE))
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='PCA Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]]))
      }
      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]])){
        withProgress(message = "Re-generating ICA plot with cluster labels", max = 1, value = 1,{
          output$plot_ica_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                useReduction = "ica",
                                showLegend = TRUE))
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='ICA Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]]))
      }

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]])){
        withProgress(message = "Re-generating tSNE plot with cluster labels", max = 1, value = 1,{
          output$plot_tsne_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "tsne",
                                                 showLegend = TRUE))
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='tSNE Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]]))
      }

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]])){
        withProgress(message = "Re-generating UMAP plot with cluster labels", max = 1, value = 1,{
          output$plot_umap_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                useReduction = "umap",
                                showLegend = TRUE))
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='UMAP Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]]))
      }
    }
    else{
      showNotification(paste0("'", input$reduction_clustering_method, "' reduction not found in input object"))
    }
  })

  #Update PCA/ICA message in clustering tab
  output$display_message_clustering <- renderText({
    if(input$reduction_clustering_method == "pca"){
      if(input$pca_significant_pc_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_slider," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{
      if(input$ica_significant_ic_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_slider," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Run tSNE
  observeEvent(input$run_tsne_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_tsne_method]])){
      withProgress(message = "Running tSNE", max = 1, value = 1, {
        vals$counts <- seuratRunTSNE(inSCE = vals$counts,
                                     useReduction = input$reduction_tsne_method,
                                     reducedDimName = "seuratTSNE",
                                     dims = input$pca_significant_pc_slider,
                                     perplexity = input$perplexity_tsne)
        vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE, tSNE = FALSE, UMAP = FALSE)
      })
      withProgress(message = "Plotting tSNE", max = 1, value = 1, {
        output$plot_tsne <- renderPlotly({
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                              useReduction = "tsne",
                              showLegend = FALSE))
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
      shinyjs::enable(selector = "div[value='Clustering']")

      showNotification("tSNE Complete")
    }
    else{
      showNotification(paste0("'", input$reduction_tsne_method, "' reduction not found in input object"))
    }
  })


  #Update PCA/ICA message in tSNE tab
  output$display_message_tsne <- renderText({
    if(input$reduction_tsne_method == "pca"){
      if(input$pca_significant_pc_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_slider," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{
      if(input$ica_significant_ic_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_slider," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Run UMAP
  observeEvent(input$run_umap_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_umap_method]])){
      withProgress(message = "Running UMAP", max = 1, value = 1, {
        vals$counts <- seuratRunUMAP(inSCE = vals$counts,
                                     useReduction = input$reduction_umap_method,
                                     reducedDimName = "seuratUMAP",
                                     dims = input$pca_significant_pc_slider,
                                     minDist = input$min_dist_umap,
                                     nNeighbors = input$n_neighbors_umap,
                                     spread = input$spread_umap)
        vals$counts <- .seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE, tSNE = FALSE, UMAP = FALSE)
      })
      withProgress(message = "Plotting UMAP", max = 1, value = 1, {
        output$plot_umap <- renderPlotly({
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                              useReduction = "umap",
                              showLegend = FALSE))
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
      shinyjs::enable(selector = "div[value='Clustering']")
      showNotification("UMAP Complete")
    }
    else{
      showNotification(paste0("'", input$reduction_umap_method, "' reduction not found in input object"))
    }
  })

  #Update PCA/ICA message in UMAP tab
  output$display_message_umap <- renderText({
    if(input$reduction_umap_method == "pca"){
      if(input$pca_significant_pc_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_slider," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{ #ICA to do
      if(input$ica_significant_ic_slider){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_slider," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Update pca significant slider maximum value with total number of computed principal components
  observe({
    req(vals$counts)
    if (!is.null(vals$counts@metadata$seurat$count_pc)) {
      updateSliderInput(session = session, inputId = "pca_significant_pc_slider", max = vals$counts@metadata$seurat$count_pc)
    }
  })

  #Update ica significant slider maximum value with total number of computed independent components
  observe({
    req(vals$counts)
    if (!is.null(vals$counts@metadata$seurat$count_ic)) {
      updateSliderInput(session = session, inputId = "ica_significant_ic_slider", max = vals$counts@metadata$seurat$count_ic)
    }
  })

  #Update tsne, umap and clustering selected number of principal components input
  observe({
    if (input$reduction_umap_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_umap_count", value = input$pca_significant_pc_slider)
    }
    else if (input$reduction_umap_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_umap_count", value = vals$counts@metadata$seurat$count_ic)
    }
    if (input$reduction_clustering_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_clustering_count", value = input$pca_significant_pc_slider)
    }
    else if (input$reduction_clustering_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_clustering_count", value = vals$counts@metadata$seurat$count_ic)
    }
    if (input$reduction_tsne_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_tsne_count", value = input$pca_significant_pc_slider)
    }
    else if (input$reduction_tsne_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_tsne_count", value = vals$counts@metadata$seurat$count_ic)
    }
  })

  #Customize heatmap (pca) with selected options
  observeEvent(input$plot_heatmap_pca_button, {
    if (!is.null(input$picker_dimheatmap_components_pca)) {
      output$plot_heatmap_pca <- renderPlot({
        seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                          dims = length(input$picker_dimheatmap_components_pca),
                          ncol = input$slider_dimheatmap_pca,
                          labels = input$picker_dimheatmap_components_pca)
      })
    }
  })

  #Customize heatmap (ica) with selected options
  observeEvent(input$plot_heatmap_ica_button, {
    if (!is.null(input$picker_dimheatmap_components_ica)) {
      output$plot_heatmap_ica <- renderPlot({
        seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                          dims = length(input$picker_dimheatmap_components_ica),
                          ncol = input$slider_dimheatmap_ica,
                          labels = input$picker_dimheatmap_components_ica)
      })
    }
  })


  #Disable tabs (THIS PART NEEDS REFACTORING)
  observe({
    if(!is.null(vals$counts)){
      shinyjs::enable(
        selector = "div[value='Normalize Data']")
      if(!is.null(vals$counts@metadata)){
        if(!is.null(vals$counts@metadata$seurat)){
          if(!"seuratScaledData" %in% assayNames(vals$counts)){
            updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "primary"))
            updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "primary"))
            shinyjs::disable(selector = "div[value='Highly Variable Genes']")
          }
          if(length(slot(vals$counts@metadata$seurat$obj, "assays")[["RNA"]]@var.features) <= 0){
            updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "primary"))
            updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "primary"))
            shinyjs::disable(selector = "div[value='Dimensionality Reduction']")
          }
          if("reductions" %in% slotNames(vals$counts@metadata$seurat$obj)){
            if(is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]])
               && is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]])){
              updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "primary"))
              updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "primary"))
              shinyjs::disable(selector = "div[value='tSNE/UMAP']")
            }
            if(is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]])
               && is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]])){
              updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "primary"))
              updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "primary"))
              shinyjs::disable(selector = "div[value='Clustering']")
            }
            if(!"seurat_clusters" %in% names(vals$counts@metadata$seurat$obj@meta.data)){
              updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "primary"))
            }
          }

        }

      }
    }
    else{
      shinyjs::disable(
        selector = "div[value='Normalize Data']")
      shinyjs::disable(
        selector = "div[value='Scale Data']")
      shinyjs::disable(
        selector = "div[value='Highly Variable Genes']")
      shinyjs::disable(
        selector = "div[value='Dimensionality Reduction']")
      shinyjs::disable(
        selector = "div[value='tSNE/UMAP']")
      shinyjs::disable(
        selector = "div[value='Clustering']")
      shinyjs::disable(
        selector = "div[value='Scale Data']")

      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='PCA Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='Elbow Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='JackStraw Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='Heatmap Plot']")

      shinyjs::disable(
        selector = ".seurat_ica_plots a[data-value='ICA Plot']")
      shinyjs::disable(
        selector = ".seurat_ica_plots a[data-value='Heatmap Plot']")

      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='PCA Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='ICA Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='tSNE Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='UMAP Plot']")
    }
  })
})
