FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV TRIMMOMATIC_VERSION=0.36
ENV PACKAGES trimmomatic=${TRIMMOMATIC_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
