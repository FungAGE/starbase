#' wiki UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @import dataspice
#'
#' @importFrom shiny NS tagList 

create_spice()
prep_attributes(data_path = file.path("data","csv"))
# edit_attributes()
# dataspice::edit_biblio()
write_spice()



mod_wiki_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    titlePanel("Metadata Wiki"),
    tabsetPanel(
      tabPanel("Metadata", DT::dataTableOutput("metadata_table"), actionButton("save", "Save")),
      tabPanel("Data", DT::dataTableOutput("data_table"))
    )
  )}
    
#' wiki Server Functions
#'
#' @noRd 
mod_wiki_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    metadata <- reactiveVal(dataspice::edit_attributes())
    
    output$metadata_table <- DT::renderDataTable({
      dat <- metadata()
      dat
    })
    
    observeEvent(input$save, {
      write.csv(metadata(), "metadata.csv")
    })
    
    output$data_table <- DT::renderDataTable(data)
  })
}