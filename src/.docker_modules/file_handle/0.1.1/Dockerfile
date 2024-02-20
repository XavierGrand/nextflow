FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV FILE_HANDLE_VERSION 0.1.1
ENV PACKAGES git=1:2.17.0* \
    ca-certificates=20180409


RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/LBMC/file_handle.git && \
    cd file_handle && \
    git checkout tags/v${FILE_HANDLE_VERSION} && \
    cp src/file_handle.py /usr/bin/file_handle.py && \
    chmod +x /usr/bin/file_handle.py
