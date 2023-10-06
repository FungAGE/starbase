#' explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_explore_ui <- function(id){
  ns <- NS(id)
            
  fluidRow(div(
    box(title="Captain Phylogeny",mod_phylogeny_ui("phylogeny_1")),
    box(title="Represented Species",mod_metadata_ui("metadata_1"))
  ))
  fluidRow(box(title="Families",
               read_tsv("/home/adrian/Systematics/Starship_Database/Starships/family/family-names.tsv") %>%
                 select(-notes) %>%
                 datatable(options = list(), class = "display",
                           callback = JS("return table;"), #rownames, colnames, container,
                           caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
                           style = "auto", width = NULL, height = NULL, elementId = NULL,
                           fillContainer = getOption("DT.fillContainer", NULL),
                           autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
                           selection = c("multiple", "single", "none"), extensions = list(),
                           plugins = NULL, editable = FALSE)
  ))
}
    
#' explore Server Functions
#'
#' @noRd 
mod_explore_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    mod_phylogeny_server("phylogeny_1")
    # mod_metadata_server("metadata_1")
  })
}
