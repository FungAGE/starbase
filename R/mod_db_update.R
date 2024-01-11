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

pool <- pool::dbPool(RSQLite::SQLite(), dbname = "SQL/starbase.sqlite")
con <- pool::poolCheckout(pool)

onStop(function() {
  pool::poolReturn(con)
  dbDisconnect(con)
})


mod_db_update_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    box(title='Update starbase entries',width=NULL,
      DTedit::dtedit_ui(ns("sql_table")))
    # shinyjs::useShinyjs(),
    # appCSS <- ".mandatory_star { color: red; }"
    # shinyjs::inlineCSS(appCSS),
    # fluidRow(
    #   actionButton(ns("add_button"), "Add", icon("plus")),
    #   actionButton(ns("edit_button"), "Edit", icon("edit")),
    #   actionButton(ns("copy_button"), "Copy", icon("copy")),
    #   actionButton(ns("delete_button"), "Delete", icon("trash-alt"))
    # ),
    # br(),
    # fluidRow(width="100%",
    #   DT::dataTableOutput(ns("sql_table"))
    # )
  )
}

#' db_update Server Functions
#'
#' @noRd
mod_db_update_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # better if object is loaded in the server
    # so make a copy of joined_ships and make reactive to inputs  
    ship_df <- joined_ships
    
    ##### Callback functions.
    appendData <- function(data, row) {
      ship_df <- rbind(data, ship_df)
      return(ship_df)
      query <- sqlAppendTable(con, "joined_ships", data, row.names = FALSE)
      dbExecute(con, query)
    }

    updateData <- function(data, olddata, row) {
      ship_df[row,] <- data[1,]
      return(ship_df)

      sql_names <- colnames(data)
      # Generate the SET part of the SQL query dynamically
      set_clause <- paste(sql_names, "= ?", collapse = ", ")

      rowid <- dbGetQuery(con, sprintf('SELECT "rowid" FROM "joined_ships" WHERE "rowid" = %d', row))

      # Construct and execute the SQL update statement
      sql_statement <- sprintf("INSERT OR REPLACE INTO %s(row_id, %s) VALUES (?, %s)", "joined_ships", set_clause, set_clause)
      params <- c(as.list(unname(data[1,])), rowid)
      dbExecute(con, sql_statement, params)
    }

    deleteData <- function(data, row) {
      ship_df <- ship_df[-row,]
      return(ship_df)
      SQL_df <- dbReadTable(con, "joined_ships")
      row_selection <- SQL_df[row, "rowid"] 
      dbExecute(con, sprintf('DELETE FROM "joined_ships" WHERE "rowid" == ("%s")', row))
    }
    
    ##### Create the DTedit object
    reactive({
      req(credentials()$user_auth)
      # TODO: add error message is not admin
      req(user_info()$permissions == "admin")
      DTedit::dtedit_server(
          id = "sql_table",
          thedata = ship_df,
          # edit.cols = c('name', 'email', 'useR', 'notes'),
          # edit.label.cols = c('Name', 'Email Address', 'Are they an R user?', 'Additional notes'),
          # input.types = c(notes='textAreaInput'),
          # view.cols = c('name', 'email', 'useR'),
          callback.update = updateData,
          callback.insert = appendData,
          callback.delete = deleteData)
    })

    # output$sql_table <- DT::renderDataTable({
    #   table<-ship_df() %>% select(-rowid)
    #   # TODO: add names/labels to table columns?
    #   # names(table) <- ...
    #   table<-datatable(table,rownames = FALSE,options = list(searching = FALSE, lengthChange = FALSE))
    # })
    
  })
}