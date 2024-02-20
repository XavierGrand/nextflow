FROM alpine:3.13.1
MAINTAINER Lauret Modolo

ENV R_PKGS R=~4.0.3 \
  R-mathlib=~4.0.3 \
  R-dev=~4.0.3 \
  R-doc=~4.0.3

ENV R_DEPS g++ \
  libxml2-dev \
  make \
  cmake \
  linux-headers \
  cairo-dev \
  libxmu-dev \
  pango-dev \
  perl \
  tiff-dev \
  icu-dev \
  libjpeg-turbo \
  pcre-dev \
  readline-dev \
  libexecinfo-dev \
  file \
  ttf-linux-libertine \
  git \
  openssl \
  autoconf \
  automake \
  libuv \
  http-parser \
  tzdata \
  libgit2-dev \
  bash

RUN apk add --update --no-cache ${R_PKGS} ${R_DEPS}
