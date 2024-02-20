FROM python:3.9-slim
MAINTAINER Lauret Modolo

ENV DEEPTOOLS_VERSION=3.5.0
RUN apt-get update -qq \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev \
        procps
RUN pip3 install deeptools==${DEEPTOOLS_VERSION}
