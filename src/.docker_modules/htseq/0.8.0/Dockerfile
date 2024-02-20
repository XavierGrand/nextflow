
FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV HTSEQ_VERSION=0.8.0
ENV PACKAGES build-essential \
             python3-pip \
             python3-setuptools \
             python3-dev \
             python3-wheel

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install numpy==1.14.3
RUN pip3 install pysam==0.15.0
RUN pip3 install HTSeq==${HTSEQ_VERSION}
