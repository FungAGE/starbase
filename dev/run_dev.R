library(shiny)
library(reactlog)
library(golem)
library(usethis)
library(here)
library(attachment)

# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Comment this if you don't want the app to be served on a random port
options(shiny.port = httpuv::randomPort())

options(shiny.reactlog=TRUE)

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
# rm(list=ls(all.names = TRUE))

# Document and reload your package
golem::document_and_reload()

# run
# tell shiny to log all reactivity
# reactlogReset()
# reactlog_enable()

# Run the application
run_app()

# pkgload::load_all()
# golem::run_dev()

# once app has closed, display reactlog from shiny
# reactlogShow()