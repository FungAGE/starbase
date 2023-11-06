#' starfish UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import processx
#' @importFrom shiny NS tagList

# TODO: add functionality for multiple inputs?
mod_starfish_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h1("Annotate Starships with starfish")),
    fluidRow(
      box(
        icon = img(src = "img/STARFISH_LOGO.png", width = "15%"),
        p("starfish is a modular toolkit for giant mobile element annotation. Built primarily for annotating Starship elements in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal Starship or a bacterial integrative and conjugative element: a 'captain' gene with zero or more 'cargo genes downstream of its 3' end. It is particularly well suited for annotating low-copy number elements in a content independent manner.")
      ),
    ),
    fluidRow(
      box(title="Run starfish on your genome",
        column(width=6,
          textInput(ns("genus"),label = "Genus name"),
          textInput(ns("species"),label = "Species name"),
          fileInput(ns("input_fasta"), "Upload a genome file in FASTA format", accept = c(".fasta", ".fa", ".fna")),
          fileInput(ns("input_gff"), "Upload annotations for this genome in GFF3 format",accept = c(".gff",".gff3"))
          ),
        column(width=6,
          checkboxInput(ns("liftover"),label = "Perform liftover annotation?",value = FALSE),
          selectInput(ns("models"),label = "Choose which gene models to annotate:",choices = c("tyr","fre","duf3723","nlr","plp"),multiple = TRUE,selected = "tyr")
               ),
        actionButton(ns("starfish"), "Launch starfish"),
        actionButton(ns("cancel"), "Cancel Run")
      )      
    )
  )
}

#' starfish Server Functions
#'
#' @noRd
mod_starfish_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$starfish, {
      validate(
        need(input$input_fasta != "", "Provide a genome in FASTA format"),
        need(input$input_gff != "", "Provide annotations in GFF3 format"),
        need(input$species != "", "Provide a genus name"),
        need(input$genus != "", "Provide a species name")
        )
      ## show spinner
      withProgress(message = "Running Starfish...", {
      myProcess <- reactiveVal(NULL)
      myProcess <- process$new("bash bin/starfish-wrapper.sh",
                               input$genus,
                               input$species,
                               fasta_in,
                               gff_in,
                               accession_in,
                               liftover,
                               models, supervise = TRUE, stdout = "")
      myProcess$wait() # wait for the process to finish
      })
    })

    # observeEvent(input$cancel, {
    #   cat(sprintf("Closing session %s\n", session$token))
    #   session$close()
    # })
    # 
    # onStop(function() {
    #   cat(sprintf("Session %s was closed\n", session$token))
    #   if (!is.null(myProcess)) {
    #     if (myProcess$is_alive()) {
    #       myProcess$kill()
    #     }
    #   }
    # })
  })
}
