FROM debian:stretch
MAINTAINER Laurent Modolo

ENV MULTIQC_VERSION=1.7
ENV PACKAGES build-essential \
             python3-pip \
             python3-setuptools \
             python3-dev \
             python3-wheel \
             procps \
             locales

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN locale-gen en_US.UTF-8
ENV LC_ALL=en_US.utf-8
ENV LANG=en_US.utf-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8


RUN pip3 install multiqc==${MULTIQC_VERSION}

