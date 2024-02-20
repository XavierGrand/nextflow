FROM alpine:3.8
MAINTAINER Laurent Modolo

ENV MUSIC_VERSION=6613c53
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES g++ \
bash \
pcre-dev \
openssl-dev \
build-base \
bzip2-dev \
xz-dev \
git \
curl

RUN apk update && \
    apk add ${PACKAGES}
RUN curl -L -o samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    http://jaist.dl.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
 && cd samtools-${SAMTOOLS_VERSION}/ \
 && ./configure --without-curses \
 && make \
 && make install

RUN git clone https://github.com/gersteinlab/MUSIC.git && \
  cd MUSIC && \
  git checkout ${MUSIC_VERSION} && \
  make clean && \
  make  && \
  cd .. && \
  mv MUSIC/bin/MUSIC /usr/bin/ && \
  mv MUSIC/bin/generate_multimappability_signal.csh /usr/bin/ && \
  mv MUSIC/bin/run_MUSIC.csh /usr/bin/ && \
  rm -Rf MUSIC

RUN chmod +x /usr/bin/*

