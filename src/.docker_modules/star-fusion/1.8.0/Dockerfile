# work from latest LTS ubuntu release
FROM ubuntu:18.04

# set the environment variables
ENV star_version 2.7.2b
ENV star_fusion_version v1.8.0
ENV samtools_version 1.9

# run update and install necessary tools
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
    libdb-dev

# install samtools dependency
WORKDIR /usr/local/bin/
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2
RUN tar -xjf /usr/local/bin/samtools-${samtools_version}.tar.bz2 -C /usr/local/bin/
WORKDIR /usr/local/bin/samtools-${samtools_version}/
RUN  ./configure
RUN  make
RUN  make install

# install perl requirements
RUN cpanm DB_File
RUN cpanm URI::Escape
RUN cpanm Set::IntervalTree
RUN cpanm Carp::Assert
RUN cpanm JSON::XS
RUN cpanm PerlIO::gzip

# download and install star aligner dependency
WORKDIR /usr/local/bin/
RUN curl -SL https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz \
    | tar -zxvC /usr/local/bin/
WORKDIR /usr/local/bin/STAR-${star_version}/source/
RUN make STAR
RUN ln -s /usr/local/bin/STAR-${star_version}/bin/Linux_x86_64/STAR /usr/local/bin/STAR
RUN ln -s /usr/local/bin/STAR-${star_version}/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# download and install star-fusion
WORKDIR /usr/local/bin
RUN curl -SL https://github.com/STAR-Fusion/STAR-Fusion/releases/download/${star_fusion_version}/STAR-Fusion-${star_fusion_version}.FULL.tar.gz \
    | tar -zxvC /usr/local/bin/
WORKDIR /usr/local/bin/STAR-Fusion-${star_fusion_version}
RUN make
WORKDIR /usr/local/bin
ENV PATH="/usr/local/bin/STAR-Fusion-${star_fusion_version}:${PATH}"

# set default command
CMD ["STAR-Fusion --help"]