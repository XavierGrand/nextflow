FROM lbmc/samtools:1.11
MAINTAINER Laurent Modolo

ENV BWA_VERSION=0.7.17
ENV SAMBLASTER_VERSION=0.1.24

ENV PACKAGES curl \
  g++

RUN apk update && \
    apk add ${PACKAGES}

RUN curl -k -L https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 -o bwa-v${BWA_VERSION}.tar.bz2 && \
tar xjf bwa-v${BWA_VERSION}.tar.bz2 && \
cd bwa-${BWA_VERSION}/ && \
make && \
cp bwa /usr/bin && \
cd .. && \
rm -R bwa-${BWA_VERSION}/

