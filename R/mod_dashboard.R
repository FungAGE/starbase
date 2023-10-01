#' dashboard UI Function
#'
#' @description main dash board to starbase
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

library(shiny)
library(shinipsum)
library(DT)

mod_dashboard_ui <- function(id){
  ns <- NS(id)
  tagList(
    h2("A Random DT"),
    DTOutput("data_table"),
    h2("A Random Plot"),
    plotOutput("plot"),
    h2("A Random Text"),
    tableOutput("text")
  )
}
    
#' dashboard Server Functions
#'
#' @noRd 
mod_dashboard_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    function(input, output, session) {
      output$data_table <- DT::renderDT({
        random_DT(5, 5)
      })
      output$plot <- renderPlot({
        random_ggplot()
      })
      output$text <- renderText({
        random_text(nwords = 50)
      })
    }
  })
}