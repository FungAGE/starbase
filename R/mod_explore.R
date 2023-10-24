#' explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom readr read_tsv
#' @importFrom shiny NS tagList 
mod_explore_ui <- function(id){
  ns <- NS(id)
  dashboardBody(
    fluidPage(
      fluidRow(div(box(title="Starship Diversity",mod_diversity_ui("diversity_1"),width = NULL),
                   box(title="Captain Phylogeny",mod_phylogeny_ui("phylogeny_1"),width = NULL))),
      fluidRow(box(title="Represented Species",mod_metadata_ui("metadata_1"),width = NULL))
    )
  )
}
    
#' explore Server Functions
#'
#' @noRd 
mod_explore_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    mod_diversity_server("diversity_1")
    mod_phylogeny_server("phylogeny_1")
    mod_metadata_server("metadata_1")
  })
}
