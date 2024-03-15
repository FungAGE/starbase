#' blast UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import DT Biostrings msaR chorddiag ggiraph stringi htmltools htmlwidgets 
#' 
#' @importFrom readr read_csv read_file
#' @importFrom shiny NS tagList

# input logic:
# 1) test query type
# 2) look for captain genes
# 3) look for cargo genes


mod_blast_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(width=8,
        box(
          title = "Input for BLAST/HMMER Searches",
          id = "blastbox",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = NULL,
          tagList(
            textAreaInput(ns("query_text"), "Paste a sequence here in FASTA format:", value = NULL, width = "1000px", height = "200px", placeholder = "Paste any nucleotide or protein sequence here"),
            fileInput(ns("query_file"), "Upload a fasta file (10 MB maximum size)",width="400px",accept = c(".fa",".fna",".fasta",".faa")),

            checkboxInput(ns("search_ship_genes"),"Screen for genes in starship sequence?",value = FALSE),
            conditionalPanel(
              condition = "input.search_ship_genes == TRUE",
              selectInput(
                ns("gene_type"), 
                "Select a specific gene model? Default is to run search for tyr genes", 
                selected = "tyr", 
                multiple = FALSE, 
                width = "200px",
                choices = c("tyr", "fre", "nlr", "DUF3723", "plp")
              ),
              ns = ns
            ),
            div(
              style = "display:inline-block",
              selectInput(ns("eval"), "Filter by e-value:", choices = c(1, 0.001, 1e-4, 1e-5, 1e-10),multiple = FALSE,selected = 1)
            ),
            rep_br(1),
            actionButton(ns("blast"), "Start BLAST/HMMER Search")
          )
        ))),
      fluidRow(
        column(width = 8,
          uiOutput(ns("ship_ui")),
          uiOutput(ns("hmm_ui")),
          uiOutput(ns("gene_ui"))
        )
      )
      # box(title="Classification",
      #   solidHeader = TRUE,
      #   collapsible = TRUE,
      #   width = NULL,
      #   uiOutput(ns("classification_ui")))
  )
}

#' blast Server Functions
#'
#' @noRd

mod_blast_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Create reactive values
    # rv <- reactiveValues()
    
    threads=4

    # gather input and set up temp file
    tmp_fasta <- tempfile(fileext = ".fa")

    observeEvent(input$blast, {
      updateBox("blastbox", action = "toggle")
    })

    submitted <- eventReactive(input$blast, {
      query_list<-parse_fasta_input(input,tmp_fasta)
      # clean query and determine sequence type
      query_list <- clean_lines(query_list)
      return(query_list)
    })

    # select correct BLAST/HMMER program
    # force the right database to be used
    # and calls the BLAST/HMMER
    # TODO: set threshold for # of returned hits?
    blastresults <- reactive({
      req(submitted())

      # ship blast
      withProgress(message = "Performing BLAST search of Starship elements in database...", {
        ship_blast_results<-run_blast(seq_type=submitted()$query_type,blast_type="ship",tmp_fasta=tmp_fasta,input.eval=input$eval,threads=threads,stitch=FALSE)
      })

      # captain family hmmsearch
      tmp_hmmer <- tempfile(fileext = ".hmmer.txt")
      withProgress(message = "Performing HMMER search of captain gene families...", {
        hmmsearch_results<-run_hmmer(submitted()$query_type,tmp_hmmer,input$eval,tmp_fasta,threads)
      })

      results<-list(ship=ship_blast_results,hmm=hmmsearch_results)

      if (length(results[["ship"]]) < 1) {
        stop("No BLAST results")
      }

      if (length(results[["hmm"]]) < 1) {
        stop("No HMMER results")
      }

      # gene blast
      if (input$search_ship_genes == TRUE) {
        gene_list <- input$gene_type
        # combine results from multiple genes
        withProgress(message = "Performing BLAST search of Starship genes...", {
          gene_blast_results<-map(gene_list, ~{
              run_blast(seq_type=submitted()$query_type,blast_type=.x,tmp_fasta=tmp_fasta,input.eval=input$eval,threads=threads,stitch=FALSE)
          })
        })
        if (length(gene_blast_results) < 1) {
          stop("No gene BLAST results")
        }
        names(gene_blast_results) <- gene_list
        results<-c(results,gene=gene_blast_results)
      }
      return(results)
    })


    ship_tabs<-reactive({
      req(blastresults()[["ship"]])
      render_output(results=blastresults()[["ship"]],
        input=input,
        name="ship",
        output=output,
        query_type=submitted()$query_type,
        plot_chord=TRUE)
      tab_content<-list(
          DT::DTOutput(ns("table_ship")),
          msaROutput(ns("alignment_ship"), width = "100%", height = "120%"))
      if (!is.null(ns("chord_ship"))) {
        tab_content <- c(chorddiag::chorddiagOutput(ns("chord_ship"), width = "100%"), tab_content)
      }
      return(tab_content)
    })

    hmm_tabs<-reactive({
      req(blastresults()[["ship"]])
      req(blastresults()[["hmm"]])
      render_output(results=blastresults()[["hmm"]],
        input=input,
        name="hmm",
        output=output,
        query_type=submitted()$query_type,
        plot_chord=FALSE)      
      tab_content<-list(
          valueBoxOutput(ns("superfam_id")),
          ggiraph::girafeOutput(ns("superfam_tree")),
          DT::DTOutput(ns("superfam_table"))
          )
      return(tab_content)
    })

    gene_tabs<-reactive({
      req(blastresults()[["ship"]])
      req(blastresults()[["gene"]])
      render_output(results=blastresults()[["gene"]][[input$gene_type]],
        input=input,
        name=input$gene_type,
        output=output,
        query_type=submitted()$query_type,
        plot_chord=FALSE)
      tab_content<-list(
          DT::DTOutput(ns(paste0("table_",input$gene_type))),
          msaROutput(ns(paste0("alignment_",input$gene_type)), width = "100%", height = "120%"))
      if (!is.null(ns(paste0("chord_",input$gene_type)))) {
        tab_content <- c(chorddiag::chorddiagOutput(ns(paste0("chord_",input$gene_type)), width = "100%"), tab_content)
      }
      return(tab_content)
    })

    output$ship_ui <- renderUI({
      req(ship_tabs())
      box(title = "Starship BLAST Results",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = NULL,
          ship_tabs()
      )
    })

    output$hmm_ui <- renderUI({
      req(hmm_tabs())
      box(title = "Captain Superfamily",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = NULL,
          hmm_tabs()
      )
    })

    output$gene_ui <- renderUI({
      req(gene_tabs())      
      box(title="Starship Gene BLAST Results",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = NULL,
        gene_tabs()
      )
    })
  })
}