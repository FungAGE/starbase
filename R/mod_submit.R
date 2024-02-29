#' submit UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import RSQLite pool shinyjs uuid dplyr stringr
#'
#' @importFrom shiny NS tagList
#' @import dplyr glue DBI 

# TODO: add checks when entries are added to database:
# i.e. sequence already present in database

mod_submit_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    textInput(ns("uploader"),label = labelMandatory("Name of curator")),
    textInput(ns("evidence"),label=labelMandatory("What evidence exists for ship annotation?")),
    textInput(ns("genus"),label = labelMandatory("Enter genus name")),
    textInput(ns("species"),label = labelMandatory("Enter species name")),
    
    # TODO: store file and put path in SQL table
    fileInput(ns("fna"),label = labelMandatory("Upload ship sequence"),accept = c(".fa",".fna",".fasta")),
    fileInput(ns("gff3"),label = "Upload ship annotations (GFF3 format)",accept = c(".gff",".gff3")),
    textAreaInput(ns("comment"), "Comment", placeholder = "", height = 100, width = "354px"),
    helpText(labelMandatory(""), paste("Mandatory field.")),
    actionButton(ns("submit_ship"), "Submit")
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