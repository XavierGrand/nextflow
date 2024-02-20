FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BCFTOOLS_VERSION=1.7
ENV PACKAGES bcftools=${BCFTOOLS_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
