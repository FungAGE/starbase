######################################
#### CURRENT FILE: DEPLOY SCRIPT #####
######################################

# Test your app

## Run checks ----
## Check the package before sending to prod
# devtools::check()

# Deploy

## Local, CRAN or Package Manager ----
## This will build a tar.gz that can be installed locally,
## sent to CRAN, or to a package manager
# devtools::build()

## Docker ----
# Creating a new Dockerfile object
library(dockerfiler)
my_dock <- Dockerfile$new(FROM="rocker/shiny-verse:4.2.3")
my_dock$MAINTAINER("Adrian Forsythe", "adrian.e.forsythe@gmail.com")
my_dock$RUN("apt-get update && apt-get upgrade -y && apt-get install -y  libglpk-dev libgmp-dev libjq-dev libsodium-dev libmagick++-dev git ncbi-blast+ hmmer python3 python3-biopython && apt-get clean && rm -rf /var/lib/apt/lists/*")
my_dock$RUN("mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/ /home/data/")
my_dock$RUN("echo \"options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)\" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site")
# TODO: there is probably a better way to automate this, i.e. using package DESCRIPTION. functions to build from renv provided in golem or dockerfiler were not flexible enough
my_dock$RUN("Rscript -e 'remotes::install_version(\"shinydashboard\",upgrade=\"never\", version = \"0.7.2\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"shinydashboardPlus\",upgrade=\"never\", version = \"2.0.3\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"shinyjs\",upgrade=\"never\", version = \"2.1.0\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"shinyalert\",upgrade=\"never\", version = \"3.0.0\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"shinipsum\", upgrade=\"never\", version = NA)'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"pool\",upgrade=\"never\", version = \"1.0.1\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"msaR\",upgrade=\"never\", version = \"0.6.0\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"JBrowseR\",upgrade=\"never\", version = \"0.10.0\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"golem\",upgrade=\"never\", version = \"0.4.1\")'")
my_dock$RUN("Rscript -e 'remotes::install_version(\"DT\",upgrade=\"never\", version = \"0.31\")'")
my_dock$RUN("Rscript -e 'remotes::install_github(\"mattflor/chorddiag\")'")
my_dock$RUN("Rscript -e 'BiocManager::install(c(\"GenomeInfoDb\",\"BiocGenerics\",\"zlibbioc\",\"S4Vectors\",\"IRanges\",\"XVector\",\"Biostrings\"),ask=F)'")

# my_dock$RUN("git clone https://github.com/FungAGE/Starships.git")

my_dock$RUN("rm -rf /srv/shiny-server/*")
my_dock$COPY(".", "/srv/shiny-server/")

my_dock$USER("shiny")
my_dock$EXPOSE(3838)
my_dock$CMD("[\"/usr/bin/shiny-server\"]")
my_dock$write()

# If you want to deploy to ShinyProxy
# golem::add_dockerfile_with_renv_shinyproxy()

## RStudio ----
## If you want to deploy on RStudio related platforms
# golem::add_rstudioconnect_file()
# golem::add_shinyserver_file()
# golem::add_shinyappsio_file()

# have to install your own app locally first
# remotes::install_local("starbase_0.0.0.9000.tar.gz")