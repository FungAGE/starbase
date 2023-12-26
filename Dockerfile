FROM rocker/shiny-verse:4.2.3
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  libglpk-dev libgmp-dev libjq-dev libsodium-dev libmagick++-dev git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("XML",upgrade="never", version = "3.99-0.16")'
RUN Rscript -e 'remotes::install_version("hunspell",upgrade="never", version = "3.0.3")'
RUN Rscript -e 'remotes::install_version("anytime",upgrade="never", version = "0.3.9")'
RUN Rscript -e 'remotes::install_version("waiter",upgrade="never", version = "0.2.5")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("crosstalk",upgrade="never", version = "1.2.1")'
RUN Rscript -e 'remotes::install_version("taxize",upgrade="never", version = "0.9.100")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.2")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("tidytree",upgrade="never", version = NA)'
RUN Rscript -e 'remotes::install_version("taxa",upgrade="never", version = "0.4.2")'
RUN Rscript -e 'remotes::install_version("spelling",upgrade="never", version = "2.2.1")'
RUN Rscript -e 'remotes::install_version("sodium",upgrade="never", version = "1.3.1")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.8.0")'
RUN Rscript -e 'remotes::install_version("shinymeta",upgrade="never", version = "0.2.0.3")'
RUN Rscript -e 'remotes::install_version("shinydashboardPlus",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("shinyalert",upgrade="never", version = "3.0.0")'
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
RUN Rscript -e 'BiocManager::install(c("GenomeInfoDb","BiocGenerics","zlibbioc","S4Vectors","IRanges","XVector","Biostrings","treeio"),ask=F)'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone /srv/shiny-server/*
COPY . /srv/shiny-server/
USER shiny
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
