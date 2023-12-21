#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard DT shinyjs sodium
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  credentials <- data.frame(
    username_id = c("myuser", "myuser1"),
    passod = sapply(c("mypass", "mypass1"), password_store),
    permission = c("Basic User", "Admin"),
    image = app_sys("img/favicon.ico"),
    stringsAsFactors = F
  )

  login <- FALSE
  USER <- reactiveValues(login = login)

  loginpage <- div(
    id = "loginpage", style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
    wellPanel(
      tags$h2("LOG IN", class = "text-center", style = "padding-top: 0;color:#333; font-weight:600;"),
      textInput("userName", placeholder = "Username", label = tagList(icon("user"), "Username")),
      passwordInput("passwd", placeholder = "Password", label = tagList(icon("unlock-alt"), "Password")),
      br(),
      div(
        style = "text-align: center;",
        actionButton("login", "SIGN IN", style = "color: white; background-color:#3c8dbc;
                                 padding: 10px 15px; width: 150px; cursor: pointer;
                                 font-size: 18px; font-weight: 600;"),
        shinyjs::hidden(
          div(
            id = "nomatch",
            tags$p("Oops! Incorrect username or password!",
              style = "color: red; font-weight: 600; 
                                            padding-top: 5px;font-size:16px;",
              class = "text-center"
            )
          )
        ),
        br(),
        br(),
        tags$code("Username: myuser  Password: mypass"),
        br(),
        tags$code("Username: myuser1  Password: mypass1")
      )
    )
  )

  observe({
    if (USER$login == FALSE) {
      if (!is.null(input$login)) {
        if (input$login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if (length(which(credentials$username_id == Username)) == 1) {
            pasmatch <- credentials["passod"][which(credentials$username_id == Username), ]
            pasverify <- password_verify(pasmatch, Password)
            if (pasverify) {
              USER$login <- TRUE
            } else {
              shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
            }
          } else {
            shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
            shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
          }
        }
      }
    }
  })

  output$logoutbtn <- renderUI({
    if (USER$login == TRUE) {
      menuItem("logout", tabName = "Logout",href = "javascript:window.location.reload(true)", icon = icon("key"))
    } else {}
  })

  output$usersidebarpanel <- renderUI({
    if (USER$login == TRUE) {
        menuItem("Update starbase Entries", tabName = "db_update", icon = icon("th"))
    } else {}
  })

  output$userloginpage<-renderUI({
    if (USER$login == TRUE) {
      Username <- isolate(input$userName)
      userBox(
        title = userDescription(
          title = credentials[which(credentials$username_id==Username),"username_id"],
          subtitle = credentials[which(credentials$username_id==Username),"permission"],
          type = 2,
          image = shinipsum::random_image(),
          # image = credentials[which(credentials$username_id==Username),"image"],
          shinipsum::random_text(nchars = 20)
        ),
        status = "warning"
        )
        } else {
          # TODO: link to login page here
          loginpage
        }
      
  })


  mod_home_server("home_1")
  mod_wiki_server("wiki_1")
  mod_explore_server("explore_1")
  mod_blast_server("blast_1")
  # mod_genome_browser_server("genome_browser_1")
  mod_blast_syn_viz_server("blast_syn_viz_1")
  mod_starfish_server("starfish_1")
  # mod_submit_server("submit_1")
  mod_db_update_server("db_update_1")
}
