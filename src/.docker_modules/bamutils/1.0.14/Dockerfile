FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BAMUTILS_VERSION=1.0.14
ENV PACKAGES git \
             build-essential \
             g++ \
             libssl-dev \
             zlib1g-dev

RUN apt-get update && apt-get -y install  ${PACKAGES}

RUN git clone https://github.com/statgen/libStatGen && \
cd libStatGen && \
git checkout fae4fca874b3b78bf9b61c0 && \
make && \
cd ../ && \
git clone https://github.com/statgen/bamUtil && \
cd bamUtil && \
git checkout v${BAMUTILS_VERSION} && \
make && \
make install
