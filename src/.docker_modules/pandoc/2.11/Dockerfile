FROM alpine:3.13
MAINTAINER Laurent Modolo

ENV PANDOC_VERSION=2.11

RUN echo "https://dl-cdn.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories \
  && apk update \
  && apk add pandoc~=${PANDOC_VERSION} make
