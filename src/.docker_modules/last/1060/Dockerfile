FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV LAST_VERSION=1060
ENV PACKAGES curl=7.58.0* \
    unzip \
    make=4.1* \
    g++ \
    zlib1g-dev=1:1.2.11* \
    ca-certificates=20180409 \
    build-essential=12.4* \
    python

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L http://last.cbrc.jp/last-${LAST_VERSION}.zip -o last-${LAST_VERSION}.zip && \
unzip last-${LAST_VERSION}.zip && \
cd last-${LAST_VERSION} && \
make && \
cp src/last* /usr/bin/ && \
cp scripts/* /usr/bin/ && \
cd .. && \
rm -Rf last-${LAST_VERSION}
