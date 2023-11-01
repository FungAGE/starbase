#' blast_viz UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

# TODO: use Biostrings here and then link to modules for making alignments, circos, dot plots, and jbrowse

mod_blast_viz_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    includeHTML(app_sys("html/kablammo.html"))
  )
}


#' blast_viz Server Functions
#'
#' @noRd 
mod_blast_viz_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    }
  )
}