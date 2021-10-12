FROM rocker/shiny-verse:4.0.3

MAINTAINER David Jenkins <dfj@bu.edu>

COPY . /sctk

RUN apt-get -y update -qq \ 
  && apt-get install -y --no-install-recommends \
    libjpeg-dev libv8-dev libbz2-dev liblzma-dev libglpk-dev libmagick++-6.q16-dev \
  && R -e "devtools::install_deps('/sctk', dependencies = TRUE)" \
  && R -e "devtools::build('/sctk')" \
  && R -e "install.packages('singleCellTK_1.7.5.tar.gz', repos = NULL, type = 'source')"

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/sctk/inst/shiny', port = 3838, host = '0.0.0.0')"]
