FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV HTSEQ_VERSION=0.11.2
ENV PACKAGES build-essential=12.4* \
             python3-pip=9.0.1* \
             python3-setuptools=39.0.1* \
             python3-dev=3.6.7* \
             python3-wheel=0.30.0*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install numpy==1.14.3
RUN pip3 install pysam==0.15.0
RUN pip3 install HTSeq==${HTSEQ_VERSION}
