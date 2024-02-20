FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV FASTQC_VERSION=0.11.5
ENV PACKAGES fastqc=${FASTQC_VERSION}* \
		perl=5.26*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
