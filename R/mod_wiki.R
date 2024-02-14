#' wiki UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import dataspice shinipsum
#'
#' @importFrom shiny NS tagList
#' @importFrom tidyr unnest

# create_spice()
# prep_attributes(data_path = file.path("data","csv"))
# prep_access(data_path = file.path("data","csv"))
#
# edit_attributes()
#
# bib<-read_tsv("Starships/MTDB/20211215.pub.dedup.db",c('ome', 'genus', 'species', 'strain', 'version', 'biosample','fna', 'faa', 'gff3','taxonomy','missing1', 'missing2','genomeSource','published','assembly_acc','acquisition_date'),na=c("","NA","Unknown","unknown")) %>%
#   distinct(ome,published) %>%
#   rename("title"="ome","citation"="published")
#
# read_csv("data/metadata/biblio.csv",col_names = T) %>%
#   bind_rows(bib) %>%
#   write_csv("data/metadata/biblio.csv")
#
# edit_biblio()
# edit_access()
# edit_creators()
#

# write_spice()
# source("bin/create_dataspice_site.R")
# temp<-whisker::whisker.render(readLines(system.file("template.html5",package = "dataspice")), jsonld_to_mustache("data/metadata/dataspice.json"))
# writeLines(temp,"data/metadata/index.html")

load_metadata<-function(){
  joined_ships %>%
    filter(!is.na(starship_family)) %>%
    mutate(starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family)) %>%
    select(starshipID, starship_family, starship_navis, starship_haplotype)
}

metadata <- load_metadata() %>%
  group_by(starship_family) %>%
  summarise(named_vec = list(starshipID)) %>%
  deframe()

mod_wiki_ui <- function(id) {
  ns <- NS(id)

  # Define the metadata
  metadata <- load_metadata()

  fluidPage(
    fluidRow(
        box(
          title = "Metadata Wiki",
          closable = FALSE,
          width = NULL,
          solidHeader = FALSE,
          collapsible = FALSE,
          p("Search through metadata for Starship families"),
          selectizeInput(
            inputId = ns("family"),
            label = "Starship Family",
            choices = unique(metadata$starship_family),
            multiple = FALSE,
            selected = NULL,
            width = "30%"
          ),
    fluidRow(
    column(width=6,fluidRow(uiOutput(ns("wiki")))),
      column(width=6,
        DT::DTOutput(ns("table"))
    ))))
  )
}

#' wiki Server Functions
#'
#' @noRd
mod_wiki_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    metadata <- load_metadata()

    observe({
      if (!is.null(input$family)) {
        metadata <- metadata %>%
          filter(starship_family %in% input$family)
      }
      output$table <- DT::renderDT({
        DT::datatable(metadata)
      })
      
    output$wiki<-renderUI({
      p(random_text(nwords = 50))
    })
      })
  })
}
