FROM alpine:3.8
MAINTAINER Laurent Modolo

ENV CUTADAPT_VERSION=1.15
ENV PACKAGES bash \
             python3 \
             python3-dev \
             musl-dev \
             gcc

RUN apk update && \
    apk add ${PACKAGES}

RUN pip3 install cutadapt==${CUTADAPT_VERSION}
