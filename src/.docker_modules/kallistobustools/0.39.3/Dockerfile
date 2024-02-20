FROM python:3.8-alpine
MAINTAINER Laurent Modolo

ENV B_VERSION=0.39.3
ENV K_VERSION=0.46.1

RUN apk add  --update --no-cache bash musl-dev linux-headers g++ cmake make build-base hdf5 hdf5-dev zlib-dev autoconf bash && \
wget https://github.com/BUStools/bustools/archive/v${B_VERSION}.tar.gz && \
tar xvf v${B_VERSION}.tar.gz && \
cd bustools-${B_VERSION} && \
mkdir build && \
cd build && \
cmake .. && \
sed -i -e 's/"Common\.hpp"/"Common\.hpp"\n#include <cmath>/g' ../src/bustools_whitelist.h && \
sed -i 's/pow/std::pow/g' ../src/bustools_whitelist.cpp && \
make && \
make install && \
wget https://github.com/pachterlab/kallisto/archive/v${K_VERSION}.tar.gz && \
tar xvf v${K_VERSION}.tar.gz && \
cd kallisto-${K_VERSION} && \
mkdir build && \
cd build && \
cmake .. && \
make && \
make install && \
wget https://github.com/BUStools/getting_started/releases/download/getting_started/t2g.py && \
chmod +x t2g.py && \
mv t2g.py /usr/local/bin/ && \
rm -R kallisto* bustools* v${K_VERSION}.tar.gz v${B_VERSION}.tar.gz

CMD ["sh"]
