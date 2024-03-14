#' genome_browser UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinyjs inlineCSS
#' @importFrom shiny NS tagList htmlOutput
#' @importFrom htmltools renderTags
#' 
#' @import JBrowseR
#'
#' @noRd

mod_genome_browser_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        selectizeInput(
            inputId = ns("ship"),
            label = "Starship ID",
                      choices = ship_ids,
                      multiple = FALSE,
                      selected = NULL,
            width = "30%"
          ),
        actionButton(ns("go"), "GO!"),
        fluidRow(inlineCSS(".form-group {margin-bottom: 0;}
                                .irs-with-grid {bottom: 0px;}
                                .irs-grid {height: 13px;}
                                .irs-grid-text {height: 0px;}
                                "),
          JBrowseROutput(ns("browserOutput"))
        )
      )
    )
  )
}

#' genome_browser Server Functions
#' @importFrom JBrowseR serve_data renderJBrowseR assembly track_feature tracks default_session JBrowseR JBrowseROutput
#'
#' @import JBrowseR
#' 
#' @noRd
mod_genome_browser_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    theme <- JBrowseR::theme("#6c81c0", "#22284c")

    # starting location for browser
    mks.range.1 <- 1
    # TODO: make this the end of the element
    mks.range.2 <- 1000

    host <- getOption("shiny.host")
    port <- 5656
    tmp_fasta <- tempfile(fileext = ".fa")
    tmp_gff <- tempfile(fileext = ".gff")

    sub_df<-eventReactive(input$go, {
      req(input$ship)
      # subset SQL data
      ship_seqs %>% filter(genome_name %in% !!input$ship)
    })

    reactive({
      header <- sub_df() %>% pull(genome_name)
      fa <- sub_df() %>% pull(genome_sequence)
      gff <- tbl(con,"genome_feature") %>% filter(genome_id %in% sub_df()$id) %>% 
        mutate(phase=".",score=".") %>%
        select(contigID,source,type,start,stop,score,strand,phase,attributes)
      chr <- gff %>% pull(contigID)
      
      # need to create tmp files for fasta, gff, and fai on the fly here in order to serve for JBrowseR
      seqinr::write_fasta(sequences=fa,names=header,tmp_fasta)
      write_tsv(gff,tmp_gff)

      data_server <- serve_data(dirname(tmp_fasta), port = port)    

      if (!is.null(tmp_gff) & !is.null(tmp_fasta) & dirname(tmp_fasta) == dirname(tmp_gff)) {
        assembly <- assembly(
          paste0(host,":", port, "/", basename(tmp_fasta)),
          bgzip = TRUE
        )

        ## create configuration for a JB2 GFF FeatureTrack
        annotations_track <- track_feature(
          paste0(host,":", port, "/", basename(tmp_gff)),
          assembly
        )

        ## create the tracks array to pass to browser
        tracks <- tracks(annotations_track)
        tracks_set <- c(annotations_track)

        default_session <- default_session(
          assembly,
          tracks_set[which(!is.null(tracks_set))]
        )
        output$browserOutput <- renderJBrowseR(JBrowseR(
          "View",
          assembly = assembly,
          tracks = tracks,
          location = paste0(chr, ":", mks.range.1, "..", mks.range.2),
          defaultSession = default_session,
          theme = theme
        ))
      }
    })

      # onStop(data_server$stop_server())
  })
}