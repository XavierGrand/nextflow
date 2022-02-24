##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="Sanity"
LABEL software.version="1.0"
LABEL software.description="Sanity"
LABEL software.website=""
LABEL software.documentation="https://github.com/jmbreda/Sanity"

##### VARIABLES #####

# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.0

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes git make g++ libgomp1 \
  && git clone https://github.com/jmbreda/Sanity.git \
  && cd Sanity \
  && cd src \
  && make \
  && cp ../bin/Sanity /usr/bin \
  && cd ../../ \
  && rm -rf Sanity \
  && apt-get remove --purge --yes git make g++ \
