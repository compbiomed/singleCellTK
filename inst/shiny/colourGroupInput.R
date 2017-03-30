colourGroupInput = function(inputId) {
  ns = NS(inputId)
  
  tagList(
    uiOutput(ns("colorChoosers"))
  )
  
}

colourGroup = function(input, output, session, options="", labels = "", value = "red", ...){
  ns = session$ns
  ids = reactive(sapply(options, function (option) paste(option, "inputId", sep="_")))

  cols <- reactive({
    sapply(1:length(options), function (i) {
      if (length(options)!=0) {
        if (!is.null(input[[as.character(ids()[i])]])){
          setNames(input[[as.character(ids()[i])]], options[i])
        }
      }
    })
  })

  output$colorChoosers = renderUI({
    if (length(options)!=0){
      L = vector("list", length(options))
      for (i in 1:length(L)){
        if (is.null(input[[as.character(ids()[i])]])) {
          color = palette()[(i %% length(palette()))+1]
        } else {
          color = input[[as.character(ids()[i])]]
        }
        L[[i]] = list(do.call(colourpicker::colourInput, 
                         list(ns(ids()[i]), label = options[i], value = color, ...)))
      }
      return(L)
    }
    
  })
  
  return(cols)
}