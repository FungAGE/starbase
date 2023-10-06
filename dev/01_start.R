########################################
############## SETUP ###################
########################################

golem::create_golem(path = "./",overwrite=TRUE,package_name="starbase",with_git=TRUE)

# add dependencies
# use_package("pkg.you.want.to.add")
# use_recommended_tests()
# use_recommended_deps()

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

## Install the required dev dependencies ----
golem::install_dev_deps()

usethis::use_mit_license()  
usethis::use_readme_rmd(open = FALSE)
usethis::use_lifecycle_badge("Experimental")## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()
golem::use_recommended_deps()

## Favicon ----
golem::use_favicon( path = "hex-starbase.png")

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)

rstudioapi::navigateToFile("dev/02_dev.R")