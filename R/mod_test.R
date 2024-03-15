#' test UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_test_ui <- function(id){
  ns <- NS(id)
  tagList(
    DT::DTOutput(ns("newtable"))
  )
}
    
#' test Server Functions
#'
#' @noRd 
mod_test_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
     output$newtable <- DT::renderDT({
      cars %>%
        as.data.frame() %>%
        DT::datatable(
          options = list(), class = "display", rownames = FALSE,
          callback = JS("return table;"), # rownames, colnames, container,
          caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
          style = "auto", width = NULL, height = NULL, elementId = NULL,
          fillContainer = getOption("DT.fillContainer", NULL),
          autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
          selection = "single", extensions = list(),
          plugins = NULL, editable = FALSE)
    })
  })
}
    
## To be copied in the UI
# mod_test_ui("test_1")
    
## To be copied in the server
# mod_test_server("test_1")
