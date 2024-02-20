FROM python:3.8-slim
MAINTAINER Lauret Modolo

ENV DEEPTOOLS_VERSION=3.5.1
RUN apt-get update -qq \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev \
        libcurl4 \
        libc6 \
        procps
RUN pip3 install pysam deeptools==${DEEPTOOLS_VERSION}
