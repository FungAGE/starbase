#' blast_viz UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_blast_viz_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' blast_viz Server Functions
#'
#' @noRd 
mod_blast_viz_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    # TODO: use Biostrings here and then link to modules for making alignments, circos, dot plots, and jbrowse
  })
}
    
## To be copied in the UI
# mod_blast_viz_ui("blast_viz_1")
    
## To be copied in the server
# mod_blast_viz_server("blast_viz_1")
