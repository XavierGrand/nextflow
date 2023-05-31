FROM rocker/r-base:4.2.3

LABEL MAINTAINER "Xavier Grand <xavier.grand@ens-lyon.fr>"

RUN apt update && apt-get install -y apt-transport-https
RUN apt install -y pandoc libfontconfig1-dev libcurl4-openssl-dev \
                   libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev \
                   libharfbuzz-dev libfribidi-dev libfreetype6-dev \
                   libpng-dev libtiff5-dev libjpeg-dev && \
  R -e "install.packages(c('devtools', 'testthat', 'roxygen2', 'BiocManager', 'plotly'))" && \
  R -e "BiocManager::install('DESeq2')"

COPY *.R .

RUN Rscript install_deseq2-wrapper.R

CMD [ "bash" ]