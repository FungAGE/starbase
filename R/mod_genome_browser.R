#' genome_browser UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_genome_browser_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' genome_browser Server Functions
#'
#' @noRd 
mod_genome_browser_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_genome_browser_ui("genome_browser_1")
    
## To be copied in the server
# mod_genome_browser_server("genome_browser_1")
