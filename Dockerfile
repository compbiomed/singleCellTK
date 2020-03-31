FROM openanalytics/r-base

#Install dependencies on Ubuntu
RUN buildDeps='libpq-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libssh2-1-dev python3-pip libv8-dev' apt-get update && apt-get install -y \
	libpq-dev \
	build-essential \
	libcurl4-openssl-dev \
	libxml2-dev \
	libssl-dev \
	libssh2-1-dev \
	libv8-dev \
	python3-pip && apt-get purge -y --auto-remove $buildDeps && echo

RUN export CFLAGS="-O3 -march=nehalem" && pip3 install scrublet virtualenv

#Add singleCellTK directory and script to docker (/tmp)
RUN mkdir -p /SCTK_docker/ && mkdir /SCTK_docker/script && mkdir /SCTK_docker/modes

ADD ./install_packages.R /SCTK_docker/script
ADD ./SCTK_runQC.R /SCTK_docker/script
ADD ./modes /SCTK_docker/modes

RUN Rscript /SCTK_docker/script/install_packages.R

ENTRYPOINT ["Rscript", "/SCTK_docker/script/SCTK_runQC.R"]