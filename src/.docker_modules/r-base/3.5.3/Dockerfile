FROM rocker/r-ver:3.5.3
MAINTAINER Remi Seraphin

RUN apt update && \
apt install libssl-dev libcurl4-openssl-dev libxml2-dev zlib1g-dev git procps
RUN R -e "install.packages('tidyverse'); \
          install.packages('ggplot2')"
