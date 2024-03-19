FROM rocker/shiny-verse:4.2.3
MAINTAINER Adrian Forsythe <adrian.e.forsythe@gmail.com>
RUN apt-get update && apt-get upgrade -y && apt-get install -y  libglpk-dev libgmp-dev libjq-dev libsodium-dev libmagick++-dev git ncbi-blast+ hmmer python3 python3-biopython && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/ /home/data/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("shinydashboardPlus",upgrade="never", version = "2.0.3")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("shinyalert",upgrade="never", version = "3.0.0")'
RUN Rscript -e 'remotes::install_version("shinipsum", upgrade="never", version = NA)'
RUN Rscript -e 'remotes::install_version("pool",upgrade="never", version = "1.0.1")'
RUN Rscript -e 'remotes::install_version("msaR",upgrade="never", version = "0.6.0")'
RUN Rscript -e 'remotes::install_version("JBrowseR",upgrade="never", version = "0.10.0")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.4.1")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.31")'
RUN Rscript -e 'remotes::install_github("mattflor/chorddiag")'
RUN Rscript -e 'remotes::install_version("ggiraph",upgrade="never", version = "0.8.8")'
RUN Rscript -e 'remotes::install_version("seqinr",upgrade="never", version = "4.2.36")'
RUN Rscript -e 'remotes::install_version("markdown",upgrade="never", version = "1.12")'
RUN Rscript -e 'remotes::install_github("YuLab-SMU/ggtree")'
RUN Rscript -e 'BiocManager::install(c("GenomeInfoDb","BiocGenerics","zlibbioc","S4Vectors","IRanges","XVector","Biostrings"),ask=F)'
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/
USER shiny
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
