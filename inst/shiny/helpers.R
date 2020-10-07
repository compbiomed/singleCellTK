# All the code in this file needs to be copied to your Shiny app, and you need
# to call `withBusyIndicatorUI()` and `withBusyIndicatorServer()` in your app.
# You can also include the `appCSS` in your UI, as the example app shows.

# =============================================

# Set up a button to have an animated loading indicator and a checkmark
# for better user experience
# Need to use with the corresponding `withBusyIndicator` server function
withBusyIndicatorUI <- function(button) {
  id <- button[["attribs"]][["id"]]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      hidden(
        img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# does the same thing as above, but put the animation inside the button label
actionButtonBusy <- function(buttonId, buttonTitle) {
  tags$div(
    `data-for-btn` = buttonId,
    actionButton(
      buttonId,
      tags$div(
        buttonTitle,
        tags$span( class = "btn-loading-container", style = "float:right",
          hidden(
            img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
            icon("check", class = "btn-done-indicator")
          )
        )
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })

  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err){
    errorFunc(err, buttonId)
  })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

appCSS <- "
.btn-loading-container {
  margin-left: 10px;
  font-size: 1.2em;
}
.btn-done-indicator {
  color: green;
}
.btn-err {
  margin-top: 10px;
  color: red;
}
"

# Accordion - formatting for collapsible section in an accordion
# Use example:
#     HTML('<div class="accordion" id="myAccordion">
#       <div class="panel">'),

#         HTML(accordionSection("1","2","myAccordion")),
#           # panel content code,
#         HTML('</div>'),

#       HTML('</div>
#     </div>')
accordionSection <- function(collapseId, panelTitle, accordionId) {
  return(
    paste(
      '<button type="button" class="btn btn-default btn-block" ',
        'data-toggle="collapse" data-target="#', collapseId, '" data-parent="#', accordionId, '">',
        panelTitle,
      '</button>
      <div id="', collapseId, '" class="collapse">',
      sep = ""
    )
  )
}

# show or hide all collapses in a list
allSections <- function(action, collapseList) {
  if (action == "hide"){
    for (i in collapseList){
      shinyjs::hide(i, anim = TRUE)
    }
  }
  else if (action == "show"){
    for (i in collapseList){
      shinyjs::show(i, anim = TRUE)
    }
  }
}

withConsoleRedirect <- function(expr) {
  options(warn = 1)
  tmpSinkfileName <- tempfile()
  tmpFD <- file(tmpSinkfileName, open = "wt")
  sink(tmpFD, type="output", split = TRUE)
  sink(tmpFD, type = "message")
  
  result <- expr
  
  sink(type = "message")
  sink()
  console.out <- readChar(tmpSinkfileName, file.info(tmpSinkfileName)$size)
  unlink(tmpSinkfileName)
  if (length(console.out) > 0) {
    insertUI(paste0("#", "console"), where = "beforeEnd",
             ui = tags$p(paste0(console.out, "\n", collapse = ""))
    )
  }
  result
}

withConsoleMsgRedirect <- function(expr) {
  withCallingHandlers({
    result <- expr
  },
  message = function(m) {
    shinyjs::html(id = "console", html = m$message, add = TRUE)
  })
  result
}


#-----------#
# Gene Sets #
#-----------#
formatGeneSetList <- function(setListStr) {
  setListArr <- strsplit(setListStr, "\n")[[1]]
  setListList <- list()
  for (set in setListArr) {
    setListList[[set]] <- set
  }
  return(setListList)
}

formatGeneSetDBChoices <- function(dbIDs, dbCats) {
  splitIds = strsplit(dbIDs, " ")
  choices <- list()
  
  for (i in seq_along(splitIds)) {
    entry <- splitIds[i][[1]]
    choices[[sprintf("%s - %s", entry, dbCats[i])]] <- entry
  }
  
  return(choices)
}


#--------------#
# QC/Filtering #
#--------------#

combineQCSubPlots <- function(output, combineP, algo, sampleList, plots, plotIds, statuses) {
  if (combineP=="all") {
    output[[plotIds[[algo]]]] <- renderPlot(plots, height = 800)
  } else {
    tabsetID <- paste0(algo, "Tabs") # for the tabsetPanel within a tab
    removeUI(selector = paste0("#", plotIds[[algo]]))
    
    for (i in seq_along(sampleList)) {
      s <- sampleList[[i]]
      sID <- paste(c(algo, s, "Tab"), collapse = "")
      subPlotID <- paste(c(algo, s, "Plot"), collapse = "")
      if (is.null(statuses[[algo]][[s]])) {
        if (i == 1) {
          appendTab(tabsetID, tabPanel(s, fluidPage(id = sID, plotOutput(outputId = subPlotID))), select = TRUE)
        } else {
          appendTab(tabsetID, tabPanel(s, fluidPage(id = sID, plotOutput(outputId = subPlotID))), select = FALSE)
        }
      }
      output[[subPlotID]] <- renderPlot(plots$Sample[[s]])
    }
  }
}


arrangeQCPlots <- function(inSCE, output, algoList, sampleList, plotIDs, statuses, redDimName) {
  uniqueSampleNames <- unique(sampleList)
  combineP <- "all"
  if (length(uniqueSampleNames) > 1) {
    combineP <- "sample"
  }
  for (a in algoList) {
    if (a == "doubletCells") {
      dcPlots <- plotDoubletCellsResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                                     reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, dcPlots, plotIDs, statuses)
    } else if (a == "cxds") {
      cxPlots <- plotCxdsResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                             reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, cxPlots, plotIDs, statuses)
    } else if (a == "bcds") {
      bcPlots <- plotBcdsResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                             reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, bcPlots, plotIDs, statuses)
    } else if (a == "cxds_bcds_hybrid") {
      cxbcPlots <- plotScdsHybridResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                                     reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, cxbcPlots, plotIDs, statuses)
    } else if (a == "decontX") {
      dxPlots <- plotDecontXResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                                reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, dxPlots, plotIDs, statuses)
    } else if (a == "QCMetrics") {
      qcmPlots <- plotRunPerCellQCResults(inSCE, sample = sampleList, combinePlot = combineP, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, qcmPlots, plotIDs, statuses)
      
    } else if (a == "scrublet") {
      sPlots <- plotScrubletResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                                reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, sPlots, plotIDs, statuses)
      
    } else if (a == "doubletFinder") {
      dfPlots <- plotDoubletFinderResults(inSCE, combinePlot = combineP, sample = sampleList, 
                                                      reducedDimName = redDimName, plotLabels = "none")
      combineQCSubPlots(output, combineP, a, uniqueSampleNames, dfPlots, plotIDs, statuses)
    }
  }
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

addToColFilterParams <- function(name, criteria, id, paramsReactive) {
  threshStr <- ""
  if (is.numeric(criteria)) {
    threshStr <- sprintf("%s > %.5f", name, criteria)
  } else {
    threshArr <- list()
    for (c in criteria) {
      threshArr <- c(threshArr, sprintf("%s == '%s'", name, c))
    }
    threshStr <- paste(threshArr, collapse = " | ")
  }
  
  entry <- list(col=name, param=threshStr, id=id)
  paramsReactive$params <- c(paramsReactive$params, list(entry))
  paramsReactive$id_count <- paramsReactive$id_count + 1
}

addToRowFilterParams <- function(name, X, Y, id, paramsReactive) {
  entry <- list(row=name, X=X, Y=Y, id=id)
  paramsReactive$params <- c(paramsReactive$params, list(entry))
  paramsReactive$id_count <- paramsReactive$id_count + 1
}


formatFilteringCriteria <- function(paramsReactive) {
  criteria = list()
  for (entry in paramsReactive) {
    criteria <- c(criteria, entry$param)
  }
  return(criteria)
}


addRowFiltersToSCE <- function(inSCE, paramsReactive) {
  for (entry in paramsReactive$params) {
    rowName <- paste0(entry$row, "_filter")
    a <- assay(inSCE, entry$row)
    vec <- rowSums(a > entry$X) > entry$Y
    rowData(inSCE)[[rowName]] <- vec
    entry$param <- sprintf("%s == T", rowName)
  }
  return(inSCE)
}




