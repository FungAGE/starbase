FROM rocker/verse:4.2.3
RUN apt-get update && apt-get install -y  libcurl4-openssl-dev libfribidi-dev libglpk-dev libgmp-dev libharfbuzz-dev libicu-dev libjq-dev libpng-dev libsodium-dev libssl-dev libtiff-dev libv8-dev libxml2-dev make pandoc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("stringi",upgrade="never", version = "1.8.3")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.5.1")'
RUN Rscript -e 'remotes::install_version("htmltools",upgrade="never", version = "0.5.7")'
RUN Rscript -e 'remotes::install_version("bslib",upgrade="never", version = "0.6.1")'
RUN Rscript -e 'remotes::install_version("rmarkdown",upgrade="never", version = "2.25")'
RUN Rscript -e 'remotes::install_version("knitr",upgrade="never", version = "1.45")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.1.4")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.3.0")'
RUN Rscript -e 'remotes::install_version("vroom",upgrade="never", version = "1.6.5")'
RUN Rscript -e 'remotes::install_version("readr",upgrade="never", version = "2.1.4")'
RUN Rscript -e 'remotes::install_version("DBI",upgrade="never", version = "1.1.3")'
RUN Rscript -e 'remotes::install_version("RColorBrewer",upgrade="never", version = "1.1-3")'
RUN Rscript -e 'remotes::install_version("htmlwidgets",upgrade="never", version = "1.6.4")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.4.4")'
RUN Rscript -e 'remotes::install_version("XML",upgrade="never", version = "3.99-0.16")'
RUN Rscript -e 'remotes::install_version("httpuv",upgrade="never", version = "1.6.13")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.8.0")'
RUN Rscript -e 'remotes::install_version("hunspell",upgrade="never", version = "3.0.3")'
RUN Rscript -e 'remotes::install_version("anytime",upgrade="never", version = "0.3.9")'
RUN Rscript -e 'remotes::install_version("waiter",upgrade="never", version = "0.2.5")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("crosstalk",upgrade="never", version = "1.2.1")'
RUN Rscript -e 'remotes::install_version("taxize",upgrade="never", version = "0.9.100")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.2")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("testthat",upgrade="never", version = NA)'
RUN Rscript -e 'remotes::install_version("tidyverse",upgrade="never", version = "2.0.0")'
RUN Rscript -e 'remotes::install_version("tidytree",upgrade="never", version = NA)'
RUN Rscript -e 'remotes::install_version("taxa",upgrade="never", version = "0.4.2")'
RUN Rscript -e 'remotes::install_version("spelling",upgrade="never", version = "2.2.1")'
RUN Rscript -e 'remotes::install_version("sodium",upgrade="never", version = "1.3.1")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.8.0")'
RUN Rscript -e 'remotes::install_version("shinymeta",upgrade="never", version = "0.2.0.3")'
RUN Rscript -e 'remotes::install_version("shinydashboardPlus",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("shinyalert",upgrade="never", version = "3.0.0")'
RUN Rscript -e 'remotes::install_version("RSQLite",upgrade="never", version = "2.3.4")'
RUN Rscript -e 'remotes::install_version("reactlog",upgrade="never", version = "1.1.1")'
RUN Rscript -e 'remotes::install_version("pool",upgrade="never", version = "1.0.1")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.3")'
RUN Rscript -e 'remotes::install_version("msaR",upgrade="never", version = "0.6.0")'
RUN Rscript -e 'remotes::install_version("metacoder",upgrade="never", version = "0.3.6")'
RUN Rscript -e 'remotes::install_version("JBrowseR",upgrade="never", version = "0.10.0")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.4.1")'
RUN Rscript -e 'remotes::install_version("ggiraph",upgrade="never", version = "0.8.8")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.31")'
RUN Rscript -e 'remotes::install_version("dataspice",upgrade="never", version = NA)'
RUN Rscript -e 'remotes::install_github("yonicd/covrpage")'
RUN Rscript -e 'remotes::install_version("colourpicker",upgrade="never", version = "1.3.0")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');library(starbase);starbase::run_app()"
