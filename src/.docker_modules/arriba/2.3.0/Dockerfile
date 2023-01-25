# work from 18.04 LTS ubuntu release
FROM ubuntu:18.04
LABEL MAINTAINER="Xavier Grand" EMAIL="xavier.grand@ens-lyon.fr"

# set the environment variables
ENV star_version 2.7.10a
ENV arriba_version v2.3.0

ARG DEBIAN_FRONTEND=noninteractive

## run update and install necessary tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    vim \
    less \
    curl \
    wget \
    libnss-sss \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    cpanminus \
    libcurl4-openssl-dev \
    libdb-dev \
    r-base \
    r-bioc-genomicalignments \
    r-bioc-genomicranges \
    libxml2

# install perl requirements
RUN cpanm DB_File
RUN cpanm URI::Escape
RUN cpanm Set::IntervalTree
RUN cpanm Carp::Assert
RUN cpanm JSON::XS
RUN cpanm PerlIO::gzip

# make dir
RUN mkdir /Xnfs
RUN mkdir /scratch
RUN mkdir /scratch/Bio
RUN chmod o+w -R /scratch/Bio

# download and install star aligner dependency
WORKDIR /usr/local/bin/
RUN curl -SL https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz \
    | tar -zxvC /usr/local/bin/
WORKDIR /usr/local/bin/STAR-${star_version}/source/
RUN make STAR
RUN ln -s /usr/local/bin/STAR-${star_version}/bin/Linux_x86_64/STAR /usr/local/bin/STAR
RUN ln -s /usr/local/bin/STAR-${star_version}/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# download and install arriba
WORKDIR /usr/local/bin/
RUN curl -SL https://github.com/suhrig/arriba/releases/download/${arriba_version}/arriba_${arriba_version}.tar.gz \
    | tar -zxvC /usr/local/bin/
WORKDIR /usr/local/bin/arriba_${arriba_version}/
RUN make

ENV PATH="/usr/local/bin/arriba_${arriba_version}/:${PATH}"
RUN echo 'alias draw_fusions.R=Rscript /usr/local/bin/arriba*/draw_fusions.R' >> ~/.bashrc

WORKDIR /home/

# set default command
CMD ["arriba"]