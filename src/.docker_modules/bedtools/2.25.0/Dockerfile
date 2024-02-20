FROM ubuntu:16.04
MAINTAINER Laurent Modolo

ENV BEDTOOLS_VERSION=2.25.0
ENV PACKAGES bedtools=${BEDTOOLS_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
