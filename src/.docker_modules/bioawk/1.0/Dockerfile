FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV BIOAWK_VERSION=1.0
ENV PACKAGES git=1:2.17* \
   build-essential=12.4* \
   ca-certificates=20180409 \
   zlib1g-dev=1:1.2.11* \
   byacc

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/lh3/bioawk.git && \
  cd bioawk && \
  git checkout tags/v${BIOAWK_VERSION} && \
  make && \
  cd .. && \
  mv bioawk/bioawk /usr/bin/ && \
  rm -Rf bioawk
