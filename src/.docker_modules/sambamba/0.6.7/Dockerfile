FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV SAMBAMBA_VERSION=0.6.7
ENV PACKAGES curl=7.58.0* \
    ca-certificates=20180409 \
    build-essential=12.4* \
    zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba_v${SAMBAMBA_VERSION}_linux.tar.bz2 -o sambamba_v${SAMBAMBA_VERSION}_linux.tar.bz2 && \
tar xvjf sambamba_v${SAMBAMBA_VERSION}_linux.tar.bz2 && \
mv sambamba /usr/bin/ && \
rm -R sambamba_v${SAMBAMBA_VERSION}_linux*
