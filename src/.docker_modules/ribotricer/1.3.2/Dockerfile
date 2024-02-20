FROM ubuntu:20.04
MAINTAINER Emmanuel Labaronne

ENV RIBOTRICER_VERSION=1.3.2
ENV PACKAGES python3\
             python3-dev\
             python3-pip

RUN apt-get update && \
    apt-get install -y ${PACKAGES} && \
    apt-get clean

RUN pip3 install ribotricer==${RIBOTRICER_VERSION}
