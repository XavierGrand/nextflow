FROM ubuntu:24.04

LABEL MAINTAINER="Xavier Grand" EMAIL="xavier.grand@inserm.fr"

# set the environment variables
ENV cellranger_version 9.0.0

## run update and install necessary tools
RUN apt update -y && apt install -y \
    build-essential \
    wget

## install cellranger

WORKDIR /opt
RUN wget -O cellranger-${cellranger_version}.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-${cellranger_version}.tar.gz?Expires=1732915202&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=aEGzxRKmYSy3NAriXQAXQ2iCkfYyazdwOXL1Hbu1Fklt8YbdkAtNaJ1ucSo8gnTLwNSG93Fj2WIVbzeSDlvOzxz4sryPvz0il3dj0pmskrAh7fYt-vmxK3Pi9gzOgbDNJIgFN1RBjbOhULm~oAUwfRSQapEMSQTVCPzd95xz1BLbxS5YG6l44b9RtKI6zZ951~ibJ2R325J2ImA2xPoIhTj~HPCKt-8azlgdiTbBlrWfMR4rFEbYkqheh48O~Heux3EGjv-5~7PQu7hRYulKR4EC9OEN8vDR2LaQD4j4GxaF9gv-3XXykaM9MLBPQYnbpzKo4YRLjri982Cz-DDIZA__"
RUN tar -xzvf cellranger-${cellranger_version}.tar.gz
ENV PATH="$PATH:/opt/cellranger-${cellranger_version}"

WORKDIR /