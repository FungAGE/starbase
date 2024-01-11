#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinydashboardPlus DT shinyjs dplyr glue shinyauthr RSQLite DBI lubridate
#' @noRd 

cookie_expiry <- 7
db <- dbConnect(SQLite(), ":memory:")
dbCreateTable(db, "sessions", c(user = "TEXT", sessionid = "TEXT", login_time = "TEXT"))
load("data/user_base.rda")
app_server <- function(input, output, session) {
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

  # call the logout module with reactive trigger to hide/show
  logout_init <- shinyauthr::logoutServer(
    id = "logout",
    active = reactive(credentials()$user_auth)
  )

  observe({
    if (credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })

  user_info <- reactive({
    credentials()$info
  })

  user_data <- reactive({
    req(credentials()$user_auth)

    if (user_info()$permissions == "admin") {
      dplyr::starwars[, 1:10]
    } else if (user_info()$permissions == "standard") {
      dplyr::storms[, 1:11]
    }
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


  output$welcome <- renderText({
    req(credentials()$user_auth)

    glue("Logged in as: {user_info()$name}")
  })

  output$loginUI <- renderUI({
    req(credentials()$user_auth)
    # BUG: first item does not render automatically after login
    tabItems(
      tabItem(
        tabName = "home",
        mod_home_ui("home_1")
      ),
      tabItem(
        tabName = "wiki",
        mod_wiki_ui("wiki_1")
      ),
      tabItem(
        tabName = "blast",
        mod_blast_ui("blast_1")
      ),
      tabItem(
        tabName = "explore",
        mod_explore_ui("explore_1")
      ),
      # tabItem(
      #   tabName = "synteny",
      #   mod_blast_syn_viz_ui("blast_syn_viz_1")
      # ),
      tabItem(
        tabName = "starfish",
        mod_starfish_ui("starfish_1")
      ),
      tabItem(
        tabName = "submit",
        mod_submit_ui("submit_1")
      ),
      tabItem(
        tabName = "db_update",
        mod_db_update_ui("db_update_1")
      )
    )
  })

  mod_home_server("home_1")
  mod_wiki_server("wiki_1")
  mod_explore_server("explore_1")
  mod_blast_server("blast_1")
  # mod_blast_syn_viz_server("blast_syn_viz_1")
  mod_starfish_server("starfish_1")
  mod_submit_server("submit_1")
  mod_db_update_server("db_update_1")
}
