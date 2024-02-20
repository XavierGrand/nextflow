FROM alpine:3.12.0
MAINTAINER Lauret Modolo

ENV R_PKGS R=~4.0.0 \
  R-mathlib=~4.0.0 \
  R-dev=~4.0.0

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
  gcc \
  gfortran \
  icu-dev \
  libjpeg-turbo \
  libpng-dev \
  openblas-dev \
  pcre-dev \
  readline-dev \
  xz-dev \
  zlib-dev \
  bzip2-dev \
  curl-dev \
  libexecinfo-dev \
  git

RUN echo "http://ftp.acc.umu.se/mirror/alpinelinux.org/v3.11/main" > /etc/apk/repositories && \
  echo "http://ftp.acc.umu.se/mirror/alpinelinux.org/v3.11/community" >> /etc/apk/repositories && \
  sed -i -e 's/v[[:digit:]]\..*\//edge\//g' /etc/apk/repositories && \
  apk add --update --no-cache ${R_PKGS} ${R_DEPS}

CMD ["R", "--no-save"]
