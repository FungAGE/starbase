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

load("data/joined_ships.rda")
table_dat<-joined_ships %>%
  filter(!is.na(starship_family)) %>%
  mutate(starship_family=gsub("fam","superfam0",starship_family),
    starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family))

mod_wiki_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
            column(width=4,
        box(title="What is a Starship?",
          width=NULL,
          collapsible=TRUE,
          tags$table(style = "width: 85%",
            tags$tr(tags$td(style = "width: 100%",
                            align = "left",
                            p("Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes) [2]. They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in Paecilomyces, cheese making in Penicillium, and enable the reansfer of formaldehyde resistance in Aspergillus nidulans and Penicillium chrysogenum."))),
            tags$tr(tags$td(style = "width: 100%",
                            align = "middle",
                            img(src="img/starship-model.png",width="85%")))))),
      column(width=4,
        box(title="starbase statistics",
          width=NULL,
          collapsible=TRUE,
          shinydashboard::valueBoxOutput(ns("total_species")),
          shinydashboard::valueBoxOutput(ns("total_ships"))
        )
      )
    ),
    fluidRow(
        box(
          title = "Metadata Wiki",
          closable = FALSE,
          width = NULL,
          solidHeader = FALSE,
          collapsible = FALSE,
          h4("Search through metadata for Starship families"),
          selectizeInput(
            inputId = ns("family"),
            label = "Starship Family",
            choices = unique(table_dat$starship_family),
            multiple = FALSE,
            selected = NULL,
            width = "30%"
          ),
    fluidRow(
        DT::DTOutput(ns("table"))
    )))
  )
}

#' wiki Server Functions
#'
#' @noRd
mod_wiki_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      if (!is.null(input$family)) {
        table_dat_sub <- table_dat %>%
          filter(starship_family %in% input$family)
      }
      output$table <- DT::renderDT({
        DT::datatable(table_dat_sub)
      })
      
    })

    output$total_ships <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(
        subtitle = "Total number of Starships in starbase",
        color="green",
        value = n_distinct(metadata$starshipID),
        icon = icon("dna")
      )
    })

    n_total_species<-metadata %>% 
      filter(!is.na(genus) & !is.na(species)) %>%
      mutate(gen_spec=paste0(genus,species)) %>% distinct(gen_spec) %>% nrow()

    output$total_species<-shinydashboard::renderValueBox({
      shinydashboard::valueBox(
        subtitle = "Fungal species with Starships",
        value = n_total_species,
        icon = icon("dna")
    )})

  })
}