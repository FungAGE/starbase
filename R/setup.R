library(golem)
library(here)

# setup

golem::create_golem(path = "./",overwrite=TRUE,package_name="starbase",with_git=TRUE)
golem::install_dev_deps()

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

golem::set_golem_options()

usethis::use_mit_license()  
usethis::use_lifecycle_badge( "Experimental" )
usethis::use_git()

golem::use_recommended_tests()
golem::use_recommended_deps()

golem::use_utils_ui()
golem::use_utils_server()
golem::use_favicon( path = "hex-starbase.png")

# Creating a module skeleton
golem::add_module(name = "dashboard") 

golem::add_js_file( "script" )
golem::add_js_handler( "handlers" )
golem::add_css_file( "custom" )

# create lorem-ipsum app

# add dependencies

# run
pkgload::load_all()
golem::run_dev()