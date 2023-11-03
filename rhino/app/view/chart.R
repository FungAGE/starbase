# app/view/chart.R

box::use(
  shiny[h3, moduleServer, NS],
)

#' @export
ui <- function(id) {
  ns <- NS(id)

  h3("Chart")
}

#' @export
server <- function(id) {
  moduleServer(id, function(input, output, session) {
    print("Chart module server part works!")
  })
}