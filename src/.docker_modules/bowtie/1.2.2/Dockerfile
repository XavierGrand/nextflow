FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BOWTIE_VERSION=1.2.2
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES bowtie=${BOWTIE_VERSION}* \
    samtools=${SAMTOOLS_VERSION}*


RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
