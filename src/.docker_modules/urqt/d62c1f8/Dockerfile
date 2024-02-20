FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV URQT_VERSION=d62c1f8
ENV PACKAGES git=1:2.17* \
   build-essential=12.4* \
   ca-certificates=20180409 \
   procps \
   zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/l-modolo/UrQt.git && \
  cd UrQt && \
  git checkout ${URQT_VERSION} && \
  make && \
  cd .. && \
  mv UrQt/UrQt /usr/bin/ && \
  rm -Rf UrQt
