FROM rocker/shiny:latest

MAINTAINER Giovanni Scala "giovanni.scala@unina.it"

# system libraries of general use
RUN apt-get update && apt-get install -y libssl-dev libxml2-dev libfontconfig1-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libpq-dev libglpk-dev


# R packages

RUN R -e "install.packages(c('remotes','BiocManager'))"

RUN R -e "BiocManager::install('BiocNeighbors')"

RUN R -e "remotes::install_github('BioinfoUninaScala/MoNETA', build_vignettes=FALSE, repos=BiocManager::repositories(),dependencies=TRUE, type='source')"

RUN R -e "remotes::install_github('daattali/shinycssloaders')"


COPY Rprofile.site /usr/local/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "MoNETA::MoNETAshiny()"]
