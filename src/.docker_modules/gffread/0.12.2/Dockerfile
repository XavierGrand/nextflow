FROM alpine:3.12
MAINTAINER Laurent Modolo

ENV GFFREAD_VERSION=0.12.2
ENV PACKAGES make \
             g++ \
             bash \
             perl

RUN apk update && \
    apk add ${PACKAGES} && \
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-${GFFREAD_VERSION}.tar.gz && \
tar -xvf gffread-${GFFREAD_VERSION}.tar.gz && \
cd gffread-${GFFREAD_VERSION}/ && \
make && \
cp gffread /usr/bin/
