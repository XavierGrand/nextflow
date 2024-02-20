FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV SRATOOLKIT_VERSION=2.8.2
ENV PACKAGES sra-toolkit=${SRATOOLKIT_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
