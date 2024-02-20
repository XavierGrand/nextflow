FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV CANU_VERSION=1.6
ENV PACKAGES canu=${CANU_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean
