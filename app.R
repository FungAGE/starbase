# Launch the ShinyApp (Do not remove this comment)
# To deploy, run: rsconnect::deployApp()
# Or use the blue button on top of this file

pkgload::load_all(export_all = FALSE,helpers = FALSE,attach_testthat = FALSE)

# Set options here
# TRUE = production mode, FALSE = development mode
options(golem.app.prod = FALSE) 
# options(shiny.testmode = TRUE)
options(shiny.autoreload = TRUE)

# Comment this if you don't want the app to be served on a random port
# options(shiny.port = httpuv::randomPort())
options('shiny.port'=3838,shiny.host='127.0.0.1')

# for bioc repos
options(repos = BiocManager::repositories())

# tell shiny to log all reactivity
options(shiny.reactlog=TRUE)

# Run the application
golex::run_app()
