library(shiny)
library(reactlog)
library(golem)
library(usethis)
library(here)
library(attachment)
# setup

create_golem(path = "./",overwrite=TRUE,package_name="starbase",with_git=TRUE)

# add dependencies
# use_package("pkg.you.want.to.add")
# use_recommended_tests()
# use_recommended_deps()

install_dev_deps()

fill_desc(
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

set_golem_options()

use_mit_license()  
use_lifecycle_badge( "Experimental" )
use_git()

use_recommended_tests()
use_recommended_deps()

use_utils_ui()
use_utils_server()
use_favicon( path = "hex-starbase.png")

# Creating a module skeleton
add_module(name = "dashboard") 
add_module(name = "home")

add_module(name = "explore")
add_module(name = "metadata")
add_module(name = "phylogeny")

add_module(name = "starfish")
add_module(name = "submit")

add_module(name = "blast") 
add_module(name = "blast_viz")
add_module(name = "shoot")
add_module(name = "genome_browser")
add_module(name = "dotplot")
add_module(name = "user")
add_module(name = "sql")
add_module(name = "db_update")

add_js_file( "script" )
add_js_handler( "handlers" )
add_css_file( "custom" )

golem::add_html_template("BlasterJS")

7# create lorem-ipsum app

# This function will read all the scripts in the R/ folder and 
# try to guess required dependencies
att_from_rscripts()

# run

# tell shiny to log all reactivity
reactlogReset()
reactlog_enable()
# or
options(shiny.reactlog=TRUE)
pkgload::load_all()
golem::run_dev()

# once app has closed, display reactlog from shiny
reactlogShow()