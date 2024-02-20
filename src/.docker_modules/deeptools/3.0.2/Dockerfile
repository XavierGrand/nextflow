FROM debian:stretch
MAINTAINER Laurent Modolo

ENV DEEPTOOLS_VERSION=3.0.2
ENV PACKAGES build-essential \
             python3-pip \
             python3-setuptools \
             python3-dev \
             python3-wheel \
             zlib1g-dev \
             procps \
             libcurl4-gnutls-dev

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install deeptools==${DEEPTOOLS_VERSION}
