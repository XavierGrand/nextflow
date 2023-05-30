FROM rocker/r-base:4.2.3
LABEL MAINTAINER "Xavier Grand <xavier.grand@ens-lyon.fr>"

RUN apt update && apt-get install -y apt-transport-https
RUN apt install -y libxml2-dev libfontconfig1-dev libssl-dev \
                   libcurl4-openssl-dev libharfbuzz-dev \
                   libfribidi-dev libfreetype6-dev libpng-dev \
                   libtiff5-dev libjpeg-dev 

## copy Rscript files
COPY ./*.R .

RUN Rscript install_pkgs.R

# command to run on container start
CMD [ "bash" ]
