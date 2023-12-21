######################################
#### CURRENT FILE: DEPLOY SCRIPT #####
######################################

# Test your app

## Run checks ----
## Check the package before sending to prod
devtools::check()

# Deploy

## Local, CRAN or Package Manager ----
## This will build a tar.gz that can be installed locally,
## sent to CRAN, or to a package manager
devtools::build()

## Docker ----
# with dependances:
golem::add_dockerfile_with_renv(extra_sysreqs = c("ncbi-blast+","hmmer","python","python3-biopython"))

## RStudio ----
## If you want to deploy on RStudio related platforms
# golem::add_rstudioconnect_file()
# golem::add_shinyserver_file()
golem::add_shinyappsio_file()

# have to install your own app locally first
remotes::install_local("starbase_0.0.0.9000.tar.gz")

# Docker ----
## If you want to deploy via a generic Dockerfile
# golem::add_dockerfile_with_renv()

# If you want to deploy to ShinyProxy
# golem::add_dockerfile_with_renv_shinyproxy()

# Deploy to Posit Connect or ShinyApps.io
# In command line.
rsconnect::deployApp(
  appName = desc::desc_get_field("Package"),
  appTitle = desc::desc_get_field("Package"),
  appFiles = c(
    # Add any additional files unique to your app here.
    "SQL/",
    "html/",
    "R/",
    "inst/",
    "tmp/",
    "data/",
    "Starships/",
    "blastdb/",
    "bin/",
    "MTDB/",
    "NAMESPACE",
    "DESCRIPTION",
    "app.R"
  ),
  appId = rsconnect::deployments(".")$appID,
  lint = FALSE,
  forceUpdate = TRUE
)
