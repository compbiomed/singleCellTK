FROM rocker/shiny-verse:latest

#Install dependencies on Ubuntu
RUN buildDeps='libpq-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libssh2-1-dev python3-pip libv8-dev pandoc' apt-get update && apt-get install -y \
	libpq-dev \
	libgeos-dev \
	build-essential \
	libcurl4-openssl-dev \
	libxml2-dev \
	libssl-dev \
	libssh2-1-dev \
	libv8-dev \
	libmagick++-dev \
	libcairo2-dev \
	pandoc \
	python3-pip && apt-get purge -y --auto-remove $buildDeps && apt-get install -y curl && echo

RUN export CFLAGS="-O3 -march=nehalem" && pip3 install --upgrade pip && pip3 install numpy llvmlite scrublet virtualenv scanpy
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-cli -y

#Add singleCellTK directory and script to docker
RUN mkdir -p /SCTK_docker/ && mkdir /SCTK_docker/script && mkdir /SCTK_docker/modes 

#ADD ./install_packages.R /SCTK_docker/script
ADD ./exec/SCTK_runQC.R ./SCTK_docker/script

#Install necessary R packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "install.packages('SeuratObject')"
RUN R -e "install.packages('scran')"
RUN R -e "install.packages('Seurat')"
#RUN R -e "install.packages('shiny')"
RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('rversions')"
RUN R -e "install.packages('usethis')"
RUN R -e "install.packages('optparse', dependencies = TRUE)"
RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('kableExtra')"
RUN R -e "BiocManager::install('TENxPBMCData')"
RUN R -e "BiocManager::install('scRNAseq')"
RUN R -e "BiocManager::install('celda')"
#RUN R -e "devtools::install_github('wleepang/shiny-directory-input')"
RUN R -e "options(timeout=360000)" \
	&& R -e "devtools::install_github('compbiomed/singleCellTK@v2.12.2', force = TRUE, dependencies = TRUE)"
RUN R -e "install.packages('Matrix', version = '1.6-1')"
RUN R -e "install.packages('irlba', type = 'source')"
RUN R -e "install.packages('SeuratObject', type = 'source')"
RUN R -e "install.packages('reticulate')"
RUN R -e "Sys.setenv(RETICULATE_PYTHON = '/usr/bin/python3')"
RUN R -e "reticulate::py_config()"

ENTRYPOINT ["Rscript", "/usr/local/lib/R/site-library/singleCellTK/exec/SCTK_runQC.R"]