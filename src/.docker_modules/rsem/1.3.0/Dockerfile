FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV RSEM_VERSION=1.3.0
ENV BOWTIE2_VERSION=2.3.4.1
ENV SAMTOOLS_VERSION=1.7
ENV PACKAGES git=1:2.17* \
   build-essential=12.4* \
   ca-certificates=20180409 \
   zlib1g-dev=1:1.2.11* \
   bowtie2=${BOWTIE2_VERSION}* \
   samtools=${SAMTOOLS_VERSION}*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/deweylab/RSEM.git RSEM_${RSEM_VERSION} && \
  cd RSEM_${RSEM_VERSION} && \
  git checkout tags/v${RSEM_VERSION} && \
  make && \
  cd .. && \
  mv RSEM_${RSEM_VERSION}/rsem-* RSEM_${RSEM_VERSION}/rsem_* RSEM_${RSEM_VERSION}/convert-* RSEM_${RSEM_VERSION}/extract-* /usr/bin/ && \
  rm -Rf RSEM_${RSEM_VERSION}
