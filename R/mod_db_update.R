#' db_update UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import DT RSQLite pool shinyjs uuid dplyr DTedit
#'
#' @importFrom shiny NS tagList
#' @import dplyr glue shinyauthr RSQLite DBI lubridate 


cookie_expiry <- 7
mod_db_update_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    shinyauthr::loginUI("login",
                  title = "Please log in to access this starbase function",
                  cookie_expiry = cookie_expiry),
                  box(title='Update starbase entries',width=NULL,DTedit::dtedit_ui(ns("sql_table"))))
}

#' db_update Server Functions
#'
#' @noRd
mod_db_update_server <- function(id) {
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
        con<-load_sql("SQL/starbase.sqlite")

        # better if object is loaded in the server
        # so make a copy of joined_ships and make reactive to inputs  
        ship_df <- joined_ships
        if(user_info()$permissions != "admin") {
          shinyalert::shinyalert(title="Insufficient permissions", type = "error",closeOnClickOutside = TRUE, closeOnEsc = TRUE)
        }
      }
    })
    
    ##### DTedit Callback functions.
    dtedit_append <- function(data, row) {
      ship_df <- rbind(data, ship_df)
      return(ship_df)
    }

    dtedit_update <- function(data, row) {
      ship_df[row,] <- data[1,]
      return(ship_df)
    }

    dtedit_delete <- function(data, row) {
      ship_df <- data[-row,]
      return(ship_df)
    }

    ##### Create the DTedit object
    reactive({
      # TODO: add error message is not admin
      req(credentials()$user_auth)
      req(user_info()$permissions == "admin")
      DTedit::dtedit_server(
          id = "sql_table",
          thedata = ship_df,
          # edit.cols = c('name', 'email', 'useR', 'notes'),
          # edit.label.cols = c('Name', 'Email Address', 'Are they an R user?', 'Additional notes'),
          # input.types = c(notes='textAreaInput'),
          # view.cols = c('name', 'email', 'useR'),
          callback.update = dtedit_update,
          callback.insert = dtedit_append,
          callback.delete = dtedit_delete)
    })

    output$sql_table <- DT::renderDataTable({
      req(credentials()$user_auth)
      req(user_info()$permissions == "admin")

      table<-ship_df %>% select(-rowid)
      # TODO: add names/labels to table columns?
      # names(table) <- ...
      table<-datatable(table,rownames = FALSE,options = list(searching = FALSE, lengthChange = FALSE))
    })
  
    reactive({
      if (credentials()$user_auth) {
        onStop(stop_sql(con))
      }
    })

  })
}