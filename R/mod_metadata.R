#' metadata UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @import tidyverse DT
#' @importFrom shiny NS tagList 
mod_metadata_ui <- function(id){
  ns <- NS(id)
  # DT::DTOutput(ns("meta_table"))
  
  # a custom table container
  sketch = htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th(colspan = 11, 'Taxonomic Information'),
        th(colspan = 9, 'Starship Information')
      ),
      tr(
        # lapply(rep(c('Length', 'Width'), 2), th)
        th("ome"),th("kingdom"),th("phylum"),th("subphylum"),th("class"),th("subclass"),th("order"),th("family"),th("genus"),th("species"),th("taxid"),th("starship_class"),th("code"),th("size"),th("tsd"),th("atir"),th("spok"),th("ars"),th("other"),th("hgt")
      )
    )
  ))

  readr::read_tsv("MTDB/joined_ships.tsv") %>%
      dplyr::select(ome,kingdom,phylum,subphylum,class,subclass,order,family,genus,species,taxid,starship_class,code,size,tsd,atir,spok,ars,other,hgt) %>%
      DT::datatable(options = list(), class = "display",rownames = FALSE,container = sketch,
                callback = JS("return table;"), #rownames, colnames, container,
                caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
                style = "auto", width = NULL, height = NULL, elementId = NULL,
                fillContainer = getOption("DT.fillContainer", NULL),
                autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
                selection = c("multiple", "single", "none"), extensions = list(),
                plugins = NULL, editable = FALSE)
  
}
    
#' metadata Server Functions
#'
#' @noRd 
mod_metadata_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  #   # a custom table container
  #   sketch = htmltools::withTags(table(
  #     class = 'display',
  #     thead(
  #       tr(
  #         th(colspan = 11, 'Taxonomic Information'),
  #         th(colspan = 9, 'Starship Information')
  #       ),
  #       tr(
  #         lapply(rep(c('Length', 'Width'), 2), th)
  #       )
  #     )
  #   ))
  # 
  #   output$meta_table<-DT::renderDT({
  #     readr::read_tsv("MTDB/joined_ships.tsv") %>% 
  #       dpylr::select(ome,kingdom,phylum,subphylum,class,subclass,order,family,genus,species,taxid,starship_class,code,size,tsd,atir,spok,ars,other,hgt) %>%
  #       DT::datatable(options = list(), class = "display",rownames = FALSE,container = sketch,
  #                 callback = JS("return table;"), #rownames, colnames, container,
  #                 caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
  #                 style = "auto", width = NULL, height = NULL, elementId = NULL,
  #                 fillContainer = getOption("DT.fillContainer", NULL),
  #                 autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
  #                 selection = c("multiple", "single", "none"), extensions = list(),
  #                 plugins = NULL, editable = FALSE)     
  #   })

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
