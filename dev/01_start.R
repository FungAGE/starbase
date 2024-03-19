########################################
############## SETUP ###################
########################################

golem::create_golem(path = "./",overwrite=TRUE,package_name="starbase",with_git=TRUE)

# add dependencies
golem::use_recommended_deps()
usethis::use_package("shiny")
usethis::use_package("shinyjs")
usethis::use_package("reactlog")
usethis::use_package("shinydashboard")
usethis::use_package("shinydashboardPlus")
usethis::use_package("shinymeta")
usethis::use_package("DBI")
usethis::use_package("glue")
usethis::use_package("RMySQL")
usethis::use_package("pool")
usethis::use_package("ggplot2")
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("readr")
usethis::use_package("plotly")
usethis::use_package("ggiraph")
usethis::use_package("htmltools")
usethis::use_package("htmlwidgets")
usethis::use_package("crosstalk")
usethis::use_package("DT")
usethis::use_package("RColorBrewer")
usethis::use_package("msaR")
usethis::use_package("stringi")
usethis::use_package("stringr")
usethis::use_package("JBrowseR")
usethis::use_package("spelling")
usethis::use_package("bslib")
usethis::use_package("shinyWidgets")
usethis::use_package("shinyalert")
usethis::use_package("colourpicker")
usethis::use_package("waiter")
usethis::use_package("vroom")
usethis::use_package("seqinr")
usethis::use_package("markdown")
# usethis::use_package("sodium")
# usethis::use_package("shinyauthr")
usethis::use_package("dockerfiler")
# usethis::use_package("webr")

# add packages with github install
remotes::install_github('jbryer/DTedit')
usethis::use_package("DTedit")
remotes::install_github('mattflor/chorddiag')
usethis::use_package("chorddiag")
remotes::install_github("YuLab-SMU/ggtree")
usethis::use_package("ggtree")
remotes::install_github('yonicd/covrpage')
usethis::use_package("covrpage")


# issue with these bioconductor packages:
renv::install(c("bioc::GenomeInfoDb","bioc::BiocGenerics","bioc::zlibbioc","bioc::S4Vectors","bioc::IRanges","bioc::XVector","bioc::Biostrings","bioc::treeio"))

# track these packages
usethis::use_package("Biostrings")
usethis::use_package("DTedit")
usethis::use_package("chorddiag")
usethis::use_package("shinipsum")

## Install the required dev dependencies ----
golem::install_dev_deps()

########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application

golem::fill_desc(
  # The Name of the package containing the App 
  pkg_name = "starbase", 
  # The Title of the package containing the App 
  pkg_title = "STARBASE", 
  # The Description of the package containing the App 
  pkg_description = "PKG_DESC.", 
  # Your First Name
  author_first_name = "Adrian", 
  # Your Last Name
  author_last_name = "Forsythe", 
  # Your Email
  author_email = "adrian.e.forsythe@gmail.com", 
  # The URL of the GitHub Repo (optional) 
  repo_url = "https://github.com/FungAGE/starbase"
)   

## Set {golem} options ----
golem::set_golem_options()

usethis::use_mit_license()  
usethis::use_readme_rmd(open = FALSE)
usethis::use_badge(badge_name="Starbase status",href="https://img.shields.io/website?url=https%3A%2F%2Fstarbase.serve.scilifelab.se%2Fapp%2Fstarbase",src="")
usethis::use_lifecycle_badge("Experimental")
## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Favicon ----
golem::use_favicon( path = "inst/app/img/favicon.png")

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)

rstudioapi::navigateToFile("dev/02_dev.R")