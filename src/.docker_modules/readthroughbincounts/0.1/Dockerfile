FROM rocker/tidyverse:4.3.1
LABEL AUTHOR="Xavier Grand"
LABEL MAINTAINER="Xavier Grand <xavier.grand@inserm.fr>"
LABEL EMAIL="<xavier.grand@ens-lyon.fr>"
LABEL current_EMAIL="<xavier.grand@inserm.fr>"
LABEL build_date="2023-10-23"

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	libcurl4-openssl-dev \
        libssl-dev \
	libxml2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libmariadb-dev \
        libbz2-dev \
	procps

## copy Rscript files

COPY src/ src/
COPY LICENSE .
COPY README.md .

## Install packages
RUN Rscript src/install_packages.R 2> install_packages.log

## Langage
CMD ["bash"]
