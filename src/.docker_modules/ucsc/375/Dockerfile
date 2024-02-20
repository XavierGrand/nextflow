FROM debian:jessie
MAINTAINER Laurent Modolo

ENV PACKAGES apt-utils \
             curl \
             build-essential \
             libssl-dev \
             libpng-dev \
             uuid-dev \
             libmysqlclient-dev \
             procps \
             rsync


RUN apt-get update && \
    apt-get install -y ${PACKAGES}

ENV UCSC_VERSION=375

RUN curl -k -L http://hgdownload.soe.ucsc.edu/admin/exe/userApps.v${UCSC_VERSION}.src.tgz -o userApps.v${UCSC_VERSION}.src.tgz &&\
tar xvf userApps.v${UCSC_VERSION}.src.tgz &&\
cd userApps/ && \
make &&\
cd .. &&\
mv userApps/bin/* /usr/bin/ &&\
rm -R userApps.v${UCSC_VERSION}.src.tgz &&\
rm -R userApps
