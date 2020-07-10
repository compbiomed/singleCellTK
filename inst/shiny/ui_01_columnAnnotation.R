shinyPanelColumnAnnotation <- fluidPage(
   bsCollapse(multiple = TRUE, open = c("Options", "Annotation"),
              bsCollapsePanel("Annotation",
                              uiOutput("output_columnAnnotation_table")
              ),
              
              bsCollapsePanel("Options",
                              tabsetPanel(
                                 tabPanel("Bins",
                                          panel(
                                             uiOutput("input_select_attribute"),
                                             uiOutput("input_select_attribute_value"),
                                             textInput("input_criteria", "criteria parameter"),
                                             selectInput("input_operator", "select comparison", choices = c("=",">","<",">=","<=")),                 
                                             textInput("input_bin_name", "bin name"),
                                             actionButton("button_confirm_bin", "Confirm Bin")
                                          )
                                 ),
                                 tabPanel("Merge Columns",
                                          panel(
                                             uiOutput("input_select_attribute_merge_1"),
                                             uiOutput("input_select_attribute_merge_2"),
                                             textInput("input_select_separator_merge", "add separator between merged values", value = "_"),
                                             actionButton("button_confirm_merge","Confirm Merge")
                                          )
                                 ),
                                 tabPanel("Magic Fill",
                                          panel(
                                             uiOutput("input_select_attribute_fill_1"),
                                             uiOutput("input_select_attribute_fill_2"),
                                             uiOutput("input_select_attribute_fill_value"),
                                             textInput("input_replace_text", "new value"),
                                             actionButton("button_confirm_fill","Confirm Fill")
                                          )
                                 ),
                                 tabPanel("Clean",
                                          panel(
                                             uiOutput("input_select_attribute_clean"),
                                             selectInput("input_removal_operation", "select removal criteria", choices = c("remove alphabets", "remove digits", "remove spaces","remove symbols")),
                                             actionButton("button_confirm_clean","Confirm Clean")
                                          )
                                 ),
                                 tabPanel("Add Column",
                                          panel(
                                             textInput("input_empty_column_name", "enter new empty column name"),
                                             textInput("input_default_value_add_column", "default value to fill"),
                                             actionButton("button_confirm_empty_column_name","Create Column")
                                          #add default value
                                          )
                                 ),
                                 tabPanel("Delete Column",
                                          panel(
                                             uiOutput("input_select_attribute_delete"),
                                             actionButton("button_confirm_delete_column", "Delete")
                                          )
                                          )
                              )
                              ),
              
              bsCollapsePanel("Save/Reset Changes", style = "success",
                              
                              fluidRow(column(6,
                                              panel("Save",
                                                    h6("Changes made to the annotation must be saved before they can be used in other modules of the toolkit:"),
                                                    actionButton("button_save_colData","Save",icon = icon("save"))
                                              )
                                                    ),
                                       column(6,
                                              panel("Reset",
                                                    h6("Reset annotation to point after changes were last saved:"),
                                                    actionButton("button_restore_phenotype","Reset",icon = icon("refresh"))
                                                    )
                                              
                                              )
                                       )
                              
                              ),
              bsCollapsePanel("Upload Annotation",
                              panel(
                                 radioGroupButtons(
                                    inputId = "colEditorChoiceRadio",
                                    label = "Select source for column annotation:",
                                    choices = c("existing colData", "upload new colData"), #update this to include values
                                    selected = "existing colData"
                                 )
                              ),
                              br(),
                                 conditionalPanel(
                                    condition = sprintf("input['%s'] == 'upload new colData'", "colEditorChoiceRadio"), #update this from above
                                    panel(
                                       h6("You can either select 'replace' to completely override the current colData if there is one available or you can select 'concatenate' to merge the new colData with the existing one."),
                                       radioGroupButtons(
                                          inputId = "colEditorUploadReplaceOption",
                                          label = "Replace or Concatenate colData:",
                                          choices = c("replace", "concatenate"),
                                          selected = "replace"
                                       ),
                                       conditionalPanel(
                                          condition = "input.colEditorUploadReplaceOption == 'concatenate'",
                                          HTML("<h6><span style='color:red'> Warning:</span> Concatenation will override the columns with same names! </h6>")
                                       ),
                                       fileInput('uploadColDataFile', 'Choose file to upload',
                                                 accept = c(
                                                    'text/csv',
                                                    'text/comma-separated-values',
                                                    '.csv'
                                                 ))
                                    )
                              )
              
                              
              )
              
              )
)