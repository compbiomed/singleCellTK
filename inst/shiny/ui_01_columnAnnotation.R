shinyPanelColumnAnnotation <- fluidPage(
   
   fluidRow(
      column(3,
             fluidRow(
                column(12,
                       panel(heading = "Upload",
                             #h4("Select Column Annotation:"),
                             radioGroupButtons(
                               inputId = "colEditorChoiceRadio",
                               label = "Select source for column annotation:",
                               choices = c("existing colData", "upload new colData"), #update this to include values
                               selected = "existing colData"
                             ),
                             conditionalPanel(
                               condition = sprintf("input['%s'] == 'upload new colData'", "colEditorChoiceRadio"), #update this from above
                               fileInput('uploadColDataFile', 'Choose file to upload',
                                         accept = c(
                                            'text/csv',
                                            'text/comma-separated-values',
                                            '.csv'
                                         ))
                             ),
                             hr(),
                             #h4("Save Changes:"),
                             h6("Changes made to the annotation must be saved before they can be used in other modules of the toolkit:"),
                             actionButton("button_save_colData","Save",icon = icon("save")),
                             br(),
                             h6("Reset annotation to point after changes were last saved:"),
                             actionButton("button_restore_phenotype","Reset",icon = icon("refresh"))
                             )
                       )
             ),
             br(),
             fluidRow(
                column(12,
                       panel(heading = "Options",
                             tabsetPanel(
                                tabPanel("Bins",
                                         uiOutput("input_select_attribute"),
                                         uiOutput("input_select_attribute_value"),
                                         textInput("input_criteria", "criteria parameter"),
                                         selectInput("input_operator", "select comparison", choices = c("=",">","<",">=","<=")),                 
                                         textInput("input_bin_name", "bin name"),
                                         actionButton("button_confirm_bin", "Confirm Bin")
                                ),
                                tabPanel("Merge Columns",
                                         uiOutput("input_select_attribute_merge_1"),
                                         uiOutput("input_select_attribute_merge_2"),
                                         actionButton("button_confirm_merge","Confirm Merge")
                                ),
                                tabPanel("Magic Fill",
                                         uiOutput("input_select_attribute_fill_1"),
                                         uiOutput("input_select_attribute_fill_2"),
                                         uiOutput("input_select_attribute_fill_value"),
                                         textInput("input_replace_text", "new value"),
                                         actionButton("button_confirm_fill","Confirm Fill")
                                ),
                                tabPanel("Clean",
                                         uiOutput("input_select_attribute_clean"),
                                         selectInput("input_removal_operation", "select removal criteria", choices = c("remove alphabets", "remove digits", "remove spaces","remove symbols")),
                                         actionButton("button_confirm_clean","Confirm Clean")
                                ),
                                tabPanel("Add Column",
                                         textInput("input_empty_column_name", "enter new empty column name"),
                                         actionButton("button_confirm_empty_column_name","Create Column")
                                         #add default value
                                )
                             )
                             )
                       )
             )
            
             ),
      column(9,
            fluidRow(
               column(12,
                      panel(heading = "Column Annotation",
                            uiOutput("output_columnAnnotation_table")
                            )
                      )
            )
            )
   ),
   
   
   
 #h2("Column Annotation"),
 #textInput("input_series_id", "Enter Series Identifier", placeholder = "GSE37751"),

 #actionButton("button_series_fetch","Fetch"),
 #textOutput("output_text"),
 #DT::dataTableOutput("output_phenotype_table"),
 #h1("Phenotype Data", align = "center"),
 #h5("description", align = "center"),
 fluidRow(
   column(4,
          # tabsetPanel(
          #   tabPanel("Bins",
          #            uiOutput("input_select_attribute"),
          #            uiOutput("input_select_attribute_value"),
          #            textInput("input_criteria", "criteria parameter"),
          #            selectInput("input_operator", "select comparison", choices = c("=",">","<",">=","<=")),                 
          #            textInput("input_bin_name", "bin name"),
          #            actionButton("button_confirm_bin", "Confirm Bin")
          #   ),
          #   tabPanel("Merge Columns",
          #            uiOutput("input_select_attribute_merge_1"),
          #            uiOutput("input_select_attribute_merge_2"),
          #            actionButton("button_confirm_merge","Confirm Merge")
          #   ),
          #   tabPanel("Magic Fill",
          #            uiOutput("input_select_attribute_fill_1"),
          #            uiOutput("input_select_attribute_fill_2"),
          #            uiOutput("input_select_attribute_fill_value"),
          #            textInput("input_replace_text", "new value"),
          #            actionButton("button_confirm_fill","Confirm Fill")
          #   ),
          #   tabPanel("Clean",
          #            uiOutput("input_select_attribute_clean"),
          #            selectInput("input_removal_operation", "select removal criteria", choices = c("remove alphabets", "remove digits", "remove spaces","remove symbols")),
          #            actionButton("button_confirm_clean","Confirm Clean")
          #   ),
          #   tabPanel("Add Column",
          #            textInput("input_empty_column_name", "enter new empty column name"),
          #            actionButton("button_confirm_empty_column_name","Create Column")
          #   ),
          #   tabPanel("Restore",
          #            h5("Warning: Any changes made to the the current phenotype data will be lost!"),
          #            actionButton("button_restore_phenotype","Confirm Restore",icon = icon("refresh"))
          #   )
          # )
   ),
   column(8,
          
          # uiOutput("output_phenotype_table")
          
   )
 )
)