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
mod_explore_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    fluidPage(title = "Explore the starbase",
      
      fluidRow(
               box(
                 title = "Starship Diversity",width=NULL,
                 valueBoxOutput(ns("total_species")),valueBoxOutput(ns("total_ships")),
                 img(
                   src = "img/heat_tree.png",
                   width = "75%",
                   style = "background-color: black;"
                 )
               )),
      fluidRow(column(width=6,
               box(
                 title = "Captain Phylogeny",width=NULL,
                 readRDS(file = "data/captain-tree.RDS"))),
      column(width=6,box(
        title = "Represented Species",width=NULL,
        DT::DTOutput(ns("meta_table")))))
      )
    )
}

#' explore Server Functions
#'
#' @noRd
mod_explore_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    table_dat <- joined_ships %>%
      dplyr::select(ome, genus, species, starship_family, starship_navis, starship_haplotype, size, dr, tir, checksum)

    # a custom table container
    sketch <- htmltools::withTags(table(
      class = "display",
      thead(
        tr(
          th(colspan = 3, "Taxonomic Information"),
          th(colspan = 7, "Starship Information")
        ),
        tr(
          lapply(colnames(table_dat), th)
        )
      )
    ))

    output$meta_table <- renderDT({
      table_dat %>%
        DT::datatable(
          options = list(), class = "display", rownames = FALSE, container = sketch,
          callback = JS("return table;"), # rownames, colnames, container,
          caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
          style = "auto", width = NULL, height = NULL, elementId = NULL,
          fillContainer = getOption("DT.fillContainer", NULL),
          autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
          selection = c("multiple", "single", "none"), extensions = list(),
          plugins = NULL, editable = FALSE
        )
    })

    output$total_ships <- renderValueBox({
    valueBox(
      value = nrow(table_dat),
      subtitle = "Total number of Starships in Database",
      icon = icon("area-chart")
    )
    })
    
    output$total_species<-renderValueBox({
      valueBox(
      value = n_distinct(paste(table_dat$genus, table_dat$species)),
      subtitle = "Total number of Represented Species",
      icon = icon("area-chart")
    )})
    
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
    #     joined_ships %>%
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
