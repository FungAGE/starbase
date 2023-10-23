#' metadata UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @import tidyverse
#' @importFrom shiny NS tagList 
mod_metadata_ui <- function(id){
  ns <- NS(id)
  renderDT(read_tsv("MTDB/joined_ships.tsv"))
    # DTOutput("tbl")
}
    
#' metadata Server Functions
#'
#' @noRd 
mod_metadata_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    # mydb <- dbConnect(RSQLite::SQLite(), "SQL/starbase.sqlite")
    # output$tbl<-
    #   dbGetQuery(mydb, 'SELECT * FROM starbase') %>%
    #     distinct(ome,genus,species,strain) %>%
    #     datatable(options = list(), class = "display",
    #               callback = JS("return table;"), #rownames, colnames, container,
    #               caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
    #               style = "auto", width = NULL, height = NULL, elementId = NULL,
    #               fillContainer = getOption("DT.fillContainer", NULL),
    #               autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
    #               selection = c("multiple", "single", "none"), extensions = list(),
    #               plugins = NULL, editable = FALSE)
    # dbDisconnect(mydb)
  })
}
