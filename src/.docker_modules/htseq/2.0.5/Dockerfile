FROM ubuntu:22.04

RUN apt-get update \
    && apt-get install -y \
     procps \
     build-essential \
     python3 \
     python3-pip

RUN pip install NumPy==2.2.6
RUN pip install pysam==0.23.3
RUN pip install matplotlib==3.10.0
RUN pip install HTSeq==2.0.5
