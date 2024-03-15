#' submit UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#'
#' @importFrom shiny NS tagList
#' @import RSQLite pool shinyjs uuid stringr glue DBI 
#' 
# TODO: add checks when entries are added to database:
# i.e. sequence already present in database
cookie_expiry <- 7
mod_submit_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(width=8,
        box(title="Submission of multiple Starships to starbase",
          width=NULL,
          status="danger",
          p("Unfortunately, we can only handle submission for one Starship at a time. If you have a batch of Starships that you'd like to submit, please send the submission via email:"),
          socialButton(href="mailto:adrian.e.forsythe@gmail.com", icon = icon("envelope"))
        ),
        box(title="Submit individual Starships to starbase",width=NULL,status="primary",
          textInput(ns("uploader"),label = labelMandatory("Name of curator")),
          textInput(ns("evidence"),label=labelMandatory("How were Starships annotated? (i.e. starfish)")),
          textInput(ns("genus"),label = labelMandatory("Enter genus name")),
          textInput(ns("species"),label = labelMandatory("Enter species name")),
          h4("Coordinates of Starship in host genome:"),
          textInput(ns("host_chr"),label = labelMandatory("Host genome contig/scaffold/chromosome ID")),
          textInput(ns("ship_start"),label = labelMandatory("Start coordinate of Starship")),
          textInput(ns("ship_end"),label = labelMandatory("End coordinate for Starship")),

          # TODO: store file and put path in SQL table
          fileInput(ns("fna"),label = labelMandatory("Upload ship sequence"),accept = c(".fa",".fna",".fasta")),
          fileInput(ns("gff3"),label = "Upload gene annotations associated with Starship sequence (GFF[3] or BED format)",accept = c(".gff",".gff3",".bed")),
          textAreaInput(ns("comment"), "Any comments about the Starship features, annotations, or host genome?", placeholder = "", height = 100, width = "354px"),
          helpText(labelMandatory(""), paste("Mandatory field.")),
          actionButton(ns("submit_ship"), "Submit")
        )
      )
    )
  )
}

#' submit Server Functions
#'
#' @noRd
mod_submit_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    con<-load_sql()
  
    # List of mandatory fields for submission
    fieldsMandatory <- c("fna","genus","species","evidence","uploader")
    fieldsAll <- c("fna","gff3","genus","species","evidence","uploader","comment")
    
    # load ship_df and make reactive to inputs  
    observeEvent(input$submit_ship,{
      # save form data into data_frame format      
      # UUIDgenerate() may be useful elsewhere
      fasta_seq<-unlist(seqinr::getSequence(seqinr::read.fasta(input$fna$datapath),as.string=TRUE))
      # TODO: parse GFF with ape::read.gff?
      gff_file <- if(is.null(input$gff3)){""} else {
        input$gff3$datapath
      }
      comments <- if(is.null(input$comment)) {""} else {input$comment}
      formData<-tibble::tibble(fna = fasta_seq,
                  gff3 = gff_file,
                  genus = input$genus,
                  species = input$species,
                  evidence = input$evidence,
                  uploader = input$uploader,
                  comment = comments)

      # define which input fields are mandatory 
      mandatoryFilled <-
        vapply(fieldsMandatory,
                function(x) {
                  !is.null(input[[x]]) && input[[x]] != ""
                },
                logical(1))
      req(all(mandatoryFilled))
      append_sql(con,"submissions",formData)
    })
    reactive({
      onStop(stop_sql(con))
    })
  })
}