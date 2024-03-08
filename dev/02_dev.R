###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Amend DESCRIPTION with dependencies read from package code parsing
## install.packages('attachment') # if needed.
attachment::att_amend_desc()

# for bioc repos
options(repos = BiocManager::repositories())

## Add modules ----
## Create a module infrastructure in R/
# Creating a module skeleton
golem::add_module(name = "home")
golem::add_module(name = "explore")
golem::add_module(name = "metadata")
golem::add_module(name = "wiki")
golem::add_module(name = "starfish")
golem::add_module(name = "submit")
golem::add_module(name = "blast") 
golem::add_module(name = "genome_browser")
golem::add_module(name = "dotplot")
golem::add_module(name = "blast_syn_viz")
golem::add_module(name = "dotplot_syn")

## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct("helpers", with_test = TRUE)
golem::add_utils("helpers", with_test = TRUE)

## External resources
## Creates .js and .css files at inst/app/www
golem::add_css_file("custom")
golem::add_sass_file("custom")

# add html templates
# golem::add_html_template("BlasterJS")
# golem::add_html_template("kablammo")

## Add internal datasets ----
# TODO: incorporate database scripts into `raw_data`?
usethis::use_data_raw(name="user_base", open=rlang::is_interactive())
usethis::use_data_raw(name="captain_tree", open=rlang::is_interactive())
usethis::use_data_raw(name="joined_ships", open=rlang::is_interactive())
usethis::use_data_raw(name="ships_with_anno", open=rlang::is_interactive())
usethis::use_data_raw(name="fna_list", open=rlang::is_interactive())
usethis::use_data_raw(name="gff_list", open=rlang::is_interactive())
usethis::use_data_raw(name="taxa_list", open=rlang::is_interactive())

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("starbase")
remotes::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
##
## (You'll need GitHub there)
# usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
usethis::use_github()

# GitHub Actions
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action("test-coverage")

# Add action for PR
usethis::use_github_action("pr-commands")

# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")