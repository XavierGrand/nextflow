FROM python:3.9-buster as build
MAINTAINER Laurent Modolo

ENV ALNTOOLS_VERSION=dd96682
ENV PACKAGES git \
    ca-certificates \
    gawk \
    procps

RUN apt-get update \
    && apt-get install -y --no-install-recommends ${PACKAGES}\
    && apt-get clean \
    && git clone https://github.com/churchill-lab/alntools.git \
    && cd alntools \
    && git checkout ${ALNTOOLS_VERSION} \
    && python setup.py install \
    && cd .. \
    && rm -R alntools \
    && pip install six \
    && apt-get autoremove --purge -y git ca-certificates

CMD ["bash"] 