#' wiki UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import dataspice
#'
#' @importFrom shiny NS tagList
#' @importFrom tidyr unnest

create_spice()
prep_attributes(data_path = file.path("data","csv"))
prep_access(data_path = file.path("data","csv"))

edit_attributes()

bib<-read_tsv("MTDB/20211215.pub.dedup.db",c('ome', 'genus', 'species', 'strain', 'version', 'biosample','fna', 'faa', 'gff3','taxonomy','missing1', 'missing2','genomeSource','published','assembly_acc','acquisition_date'),na=c("","NA","Unknown","unknown")) %>%
  distinct(ome,published) %>%
  rename("title"="ome","citation"="published")

read_csv("data/metadata/biblio.csv",col_names = T) %>%
  bind_rows(bib) %>%
  write_csv("data/metadata/biblio.csv")

edit_biblio()
edit_access()
edit_creators()

write_spice()
build_site(out_path = "data/metadata/index.html")

library(jsonlite)
metadata <- fromJSON(txt = "data/metadata/dataspice.json")

metadata <- read_tsv("MTDB/joined_ships.tsv") %>%
  filter(!is.na(starship_family)) %>%
  mutate(starship_family=ifelse(grepl("^fam",starship_family) & !is.na(code),code,starship_family)) %>%
  select(starshipID,starship_family,starship_navis,starship_haplotype) %>%
  group_by(starship_family) %>%
  summarise(named_vec = list(starshipID)) %>%
  deframe()

mod_wiki_ui <- function(id) {
  ns <- NS(id)
  # Define the metadata
  metadata <- read_tsv("MTDB/joined_ships.tsv") %>%
    filter(!is.na(starship_family)) %>%
    mutate(starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family)) %>%
    select(starshipID, starship_family, starship_navis, starship_haplotype)


  fluidPage(
    fluidRow(
      column(
        3,
    box(
        title = "Metadata Wiki", 
        closable = FALSE, 
        width = 12,
        solidHeader = FALSE, 
        collapsible = FALSE,
        selectizeInput(
          inputId = ns("family"),
          label = "Starship Family",
          choices = unique(metadata$starship_family),
          multiple = TRUE,
          selected = NULL,
          width = "100%"
        )
      )
  ),
  column(9,
      DT::DTOutput(ns("table"))
    )
  )
  )
}

#' wiki Server Functions
#'
#' @noRd
mod_wiki_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Define the metadata
    metadata <- read_tsv("MTDB/joined_ships.tsv") %>%
      filter(!is.na(starship_family)) %>%
      mutate(starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family)) %>%
      select(starshipID, starship_family, starship_navis, starship_haplotype)

    ##### Create the DTedit object
    observe({
      dat <- metadata %>%
        filter(starship_family == input$family)
      output$table <- DT::renderDT({
        DT::datatable(dat)
      })
    })
  })
}
