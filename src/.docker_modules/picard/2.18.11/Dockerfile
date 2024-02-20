FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV PICARD_VERSION=2.18.11
ENV PACKAGES default-jre=2:1.1* \
    curl=7.58.0* \
    ca-certificates=20180409


RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -o picard.jar && \
  mkdir -p /usr/share/java/  && \
  mv picard.jar /usr/share/java/

COPY PicardCommandLine /usr/bin/
RUN chmod +x /usr/bin/PicardCommandLine
