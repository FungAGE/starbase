#' explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import RSQLite pool
#' @importFrom readr read_tsv
#' @importFrom shiny NS tagList

load("data/joined_ships.rda")

load_metadata<-function(){
  joined_ships %>%
    filter(!is.na(starship_family)) %>%
    mutate(starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family)) %>%
    select(starshipID, starship_family, starship_navis, starship_haplotype)
}

metadata <- load_metadata() %>%
  group_by(starship_family) %>%
  summarise(named_vec = list(starshipID)) %>%
  tibble::deframe()

mod_explore_ui <- function(id) {
  ns <- NS(id)
  fluidPage(title = "Explore the starbase",
    fluidRow(
              box(
                title = "Starship Diversity",width=12,
                valueBoxOutput(ns("total_species")),valueBoxOutput(ns("total_ships")),
                img(
                  src = "img/heat_tree.png",
                  width = "50%",
                  style = "background-color: black;"
                ),br()
              )),
    fluidRow(box(
                title = "Captain Phylogeny and Represented Species",width=NULL,
              column(width=4,
                h4("Select tyr Family: "),
                selectizeInput(
                  inputId = ns("family"),
                  label = "Starship Family",
                  choices = unique(metadata$starship_family),
                  multiple = FALSE,
                  selected = NULL,
                  width = "30%"
                ),
                ggiraph::girafeOutput(ns("captain_tree"))),
              column(width=6,
                fluidRow(valueBoxOutput(ns("total_family_ships"))),
                fluidRow(DT::DTOutput(ns("meta_table")))
              ))))}

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

    output$total_ships <- renderValueBox({
    valueBox(
      value = nrow(table_dat),
      subtitle = "Total number of Starships in starbase",
      icon = icon("dna")
    )
    })

    output$total_species<-renderValueBox({
      valueBox(
      value = n_distinct(paste(table_dat$genus, table_dat$species)),
      subtitle = "Fungi represented in starbase",
      icon = icon("dna")
    )})
    
    gg<-readRDS(file = "data/captain-tree.RDS")

    selected_state <- reactive({
      input$captain_tree_selected
    })

    output$console <- renderPrint({
      input$captain_tree_hovered
    })

    output$captain_tree <- renderGirafe({gg})

    observeEvent(input$reset, {
      session$sendCustomMessage(type = 'plot_set', message = character(0))
    })

    output$total_family_ships <- renderValueBox({
      valueBox(
        value = table_dat %>% filter(starship_family %in% gsub("-","",selected_state())) %>% nrow(),
        subtitle = paste0("Total number of ",gsub("-","",selected_state())," Starships in Database"),
        icon = icon("dna")
      )
    })

    output$meta_table <- renderDT({
      tab_dat<-table_dat %>%
        filter(starship_family %in% gsub("-","",selected_state()))
      if( nrow(tab_dat) < 1 ) return(NULL)
      
      tab_dat %>%
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

    # Using a local cache
    # cache_dir <- cache_filesystem("cache")
    # We memoise the dbGetQuery,
    #  so that every time this function is called with
    # the same parameters, 
    # the SQL query is not actually run,
    #  but the results are fetched from the cache 
    # m_get_query <- memoise(DBI::dbGetQuery, cache=cache_dir)

    # output$tbl<-
    #   m_get_query(mydb, 'SELECT * FROM joined_ships') %>%
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
