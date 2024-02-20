FROM alpine:3.8
MAINTAINER Laurent Modolo

ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES git \
             make \
             gcc \
             musl-dev \
             zlib-dev \
             ncurses-dev \
             bzip2-dev \
             xz-dev \
             bash

RUN apk update && \
    apk add ${PACKAGES}

RUN git clone https://github.com/samtools/htslib.git && \
cd htslib && \
git checkout ${SAMTOOLS_VERSION}  && \
cd .. && \
git clone https://github.com/samtools/samtools.git && \
cd samtools && \
git checkout ${SAMTOOLS_VERSION} && \
make && \
cp samtools /usr/bin/
