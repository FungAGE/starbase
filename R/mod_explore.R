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
    fluidPage(
      fluidRow(
        column(width = 6, box(title = paste("Total number of Starships in Database:", n_distinct(joined_ships$starshipID)))),
        column(width = 6, box(title = paste("Total number of Represented Species:", n_distinct(paste(joined_ships$genus, joined_ships$species)))))
      ),
      fluidRow(
        column(
          width = 6,
          box(
            title = "Starship Diversity",
            img(
              src = "img/heat_tree.png",
              width = "100%",
              style = "background-color: black;"
            )
          )
        ),
        column(width = 6, box(
          title = "Captain Phylogeny",
          readRDS(file = "data/captain-tree.RDS")
        ))
      ),
      fluidRow(box(
        title = "Represented Species",
        DT::DTOutput(ns("meta_table"))
      ))
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

    output$meta_table <- renderDataTable({
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
