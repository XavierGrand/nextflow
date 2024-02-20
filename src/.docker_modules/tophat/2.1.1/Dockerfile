FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV TOPHAT_VERSION=2.1.1
ENV PACKAGES tophat=${TOPHAT_VERSION}*\
             bowtie=1.2.2*\
             libsys-hostname-long-perl

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES}  && \
    apt-get clean
