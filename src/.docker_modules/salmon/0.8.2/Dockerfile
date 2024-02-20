FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV SALMON_VERSION=0.8.2
ENV PACKAGES curl=7.58.0* \
   ca-certificates=20180409

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L  https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz -o Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
tar xzf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
mv Salmon-${SALMON_VERSION}_linux_x86_64/bin/* /usr/bin/ && \
rm -Rf Salmon-${SALMON_VERSION}*
