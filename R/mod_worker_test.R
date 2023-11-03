#' worker_test UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

worker <- shiny.worker::initialize_worker()

mod_worker_test_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      
      # Application title
      titlePanel("shiny.worker demo"),
      
      # Sidebar with a slider input for number of bins
      sidebarLayout(
        sidebarPanel(
          div("Play with the slider. Histogram will be still responsive, even if job is running:"),
          br(),
          sliderInput(ns("bins"),
                      "Number of bins:",
                      min = 1,
                      max = 50,
                      value = 30),
          div("Then try to run new job again:"),
          br(),
          actionButton("triggerButton", "Run job (5 sec.)")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          fluidRow(
            column(6, plotOutput("distPlot")),
            column(6, plotOutput("FuturePlot")))
        )
      )
    )
  )
}
    
#' worker_test Server Functions
#'
#' @noRd 
mod_worker_test_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2]
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
    
    plotValuesPromise <- worker$run_job("plotValuesPromise", function(args) {
      Sys.sleep(5)
      cbind(rnorm(1000), rnorm(1000))
    },
    args_reactive = reactive({
      input$triggerButton
      print("triggered!")
      ""
    })
    )
    
    output$FuturePlot <- renderPlot({
      x <- plotValuesPromise()
      title <- if (is.null(x$result)) "Job is running..." else "There you go"
      points <- if (is.null(x$result)) cbind(c(0), c(0)) else x$result
      plot(points, main = title)
    })
  })
}
    
## To be copied in the UI
# 
    
## To be copied in the server
# mod_worker_test_server("worker_test_1")