FROM python:3.9-slim
MAINTAINER Laurent Modolo

ENV MACS3_VERSION=3.0.0a6
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
RUN pip install macs3==${MACS3_VERSION}
