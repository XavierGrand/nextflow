FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BOWTIE2_VERSION=2.3.4.1
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES bowtie2=${BOWTIE2_VERSION}* \
    samtools=${SAMTOOLS_VERSION}* \
    perl=5.26.1*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
