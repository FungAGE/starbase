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
#' @import dplyr glue shinyauthr RSQLite DBI lubridate 

# TODO: add checks when entries are added to database:
# i.e. sequence already present in database
cookie_expiry <- 7
mod_submit_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
  shinyauthr::loginUI("login",
                title = "Please log in to access this starbase function",
                cookie_expiry = cookie_expiry),
                uiOutput(ns("submission_ui")))
}

#' submit Server Functions
#'
#' @noRd
mod_submit_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    db <- dbConnect(SQLite(), ":memory:")
    dbCreateTable(db, "sessions", c(user = "TEXT", sessionid = "TEXT", login_time = "TEXT"))
    load("data/user_base.rda")

    # call login module supplying data frame, user and password cols and reactive trigger
    credentials <- shinyauthr::loginServer(
      id = "login",
      data = user_base,
      user_col = user,
      pwd_col = password_hash,
      sodium_hashed = TRUE,
      cookie_logins = TRUE,
      sessionid_col = sessionid,
      cookie_getter = get_sessions_from_db,
      cookie_setter = add_session_to_db,
      log_out = reactive(logout_init())
    )

    logout_init <- shinyauthr::logoutServer(
      id = "logout",
      active = reactive(credentials()$user_auth)
    )
    
    user_info <- reactive({
      credentials()$info
    })

    output$welcome <- renderText({
      req(credentials()$user_auth)

      glue("Logged in as {user_info()$permissions} user: {user_info()$name}")
    })

    # TODO: add user page
    output$user_page<-renderUI({
      req(credentials()$user_auth)

      fluidRow(
        box(
          width = 12,
          tags$h2(glue("Your permission level is: {user_info()$permissions}.
                      You logged in at: {user_info()$login_time}.")),
            width = NULL,
            status = "primary"
          )
      )
    })

    observe({
      if (credentials()$user_auth) {
        con<-load_sql("../Starships/SQL/starbase.sqlite")
        if(user_info()$permissions != "admin") {
          shinyalert::shinyalert(title="Insufficient permissions", type = "error",closeOnClickOutside = TRUE, closeOnEsc = TRUE)
        }
      }
    })
  
    # List of mandatory fields for submission
    fieldsMandatory <- c("fna","genus","species","evidence","uploader")
    fieldsAll <- c("fna","gff3","genus","species","evidence","uploader","comment")
    
    # load ship_df and make reactive to inputs  
    observeEvent(input$submit_ship,{
      req(credentials()$user_auth)
      req(user_info()$permissions == "admin")

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

    output$submission_ui<-renderUI({
      req(credentials()$user_auth)
      req(user_info()$permissions == "admin")

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
    })
    reactive({
      if (credentials()$user_auth) {
        onStop(stop_sql(con))
      }
    })
  })
}