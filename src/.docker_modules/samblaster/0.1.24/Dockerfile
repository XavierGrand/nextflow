FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV SAMBLASTER_VERSION=0.1.24
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES curl=7.58.0* \
    samtools=${SAMTOOLS_VERSION}* \
    ca-certificates=20180409 \
    build-essential=12.4* \
    zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/GregoryFaust/samblaster/releases/download/v.${SAMBLASTER_VERSION}/samblaster-v.${SAMBLASTER_VERSION}.tar.gz -o samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
tar xvf samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
cd samblaster-v.${SAMBLASTER_VERSION}/ && \
make && \
cp samblaster /usr/bin && \
cd .. && \
rm -R samblaster-v.${SAMBLASTER_VERSION}/
