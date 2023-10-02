#' dotplot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_dotplot_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' dotplot Server Functions
#'
#' @noRd 
mod_dotplot_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_dotplot_ui("dotplot_1")
    
## To be copied in the server
# mod_dotplot_server("dotplot_1")
