FROM python:3.9-slim
MAINTAINER Laurent Modolo

ENV G2GTOOLS_VERSION=0.2.8

RUN apt update \
 && apt install -y wget build-essential zlib1g-dev libbz2-dev liblzma-dev procps \
 && wget https://github.com/churchill-lab/g2gtools/archive/v${G2GTOOLS_VERSION}.tar.gz \
 && tar -xvf v${G2GTOOLS_VERSION}.tar.gz \
 && cd g2gtools-${G2GTOOLS_VERSION} \
 && pip install numpy \
 && pip install pysam \
 && make install 


CMD ["bash"]
