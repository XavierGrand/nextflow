FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV KALLISTO_VERSION=0.43.1
ENV PACKAGES curl=7.58.0* \
   ca-certificates=20180409

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L  https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VERSION}/kallisto_linux-v${KALLISTO_VERSION}.tar.gz -o kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
tar xzf kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
cp kallisto_linux-v${KALLISTO_VERSION}/kallisto /usr/bin && \
rm -Rf kallisto_linux-v${KALLISTO_VERSION}*
