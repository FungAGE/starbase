# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# options(shiny.testmode = TRUE)
options(shiny.autoreload = TRUE)

# host address
options(shiny.host = "127.0.0.1")

# address accessible by other computers on the same network
# options(shiny.host = "130.238.82.67")

# accept any connection (not just from localhost)
# options(shiny.host = "0.0.0.0")

# port
# can assume any value that you want (just assure to avoid to select ports used by other services like ssh or http)
options(shiny.port = 5858)
# or choose a random port
# options(shiny.port = httpuv::randomPort())

options(shiny.launch.browser = FALSE)
# options(quiet = TRUE)

# tell shiny to log all reactivity
options(shiny.reactlog=TRUE)

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

# Run the application
starbase::run_app()
