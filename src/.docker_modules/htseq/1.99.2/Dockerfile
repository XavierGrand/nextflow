FROM ubuntu:20.04
LABEL maintainer= Laurent Modolo

RUN apt-get update \
    && apt-get install -y \
     procps \
     build-essential \
     python3 \
     python3-pip

RUN pip install HTSeq==1.99.2