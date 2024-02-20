FROM rocker/tidyverse:4.3.1
LABEL AUTHOR="Xavier Grand"
LABEL MAINTAINER="Xavier Grand <xavier.grand@inserm.fr>"
LABEL EMAIL="<xavier.grand@inserm.fr>"
LABEL build_date="2024-02-20"

RUN apt update \
	&& apt install -y --no-install-recommends \
	libcurl4-openssl-dev \
    libssl-dev \
	libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libmariadb-dev \
    libbz2-dev \
	procps

## run update and install necessary tools
RUN apt update -y && apt install -y \
    build-essential \
    git \
    rna-star \
    samtools

RUN git clone https://github.com/PierreBSC/Viral-Track.git

COPY install_packages.R .

RUN Rscript install_packages.R

WORKDIR /Viral-Track/

CMD [ "bash" ]