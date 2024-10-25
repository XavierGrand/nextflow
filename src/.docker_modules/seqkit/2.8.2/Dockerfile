FROM ubuntu:24.04

#Begin: install prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        git \
        libcurl3-dev \
        libfreetype6-dev \
        libpng12-dev \
        libzmq3-dev \
        locate \
        pkg-config \
        rsync \
        software-properties-common \
        sudo \
        unzip \
        wget \
        zip \
        zlib1g-dev \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
#End: install prerequisites

#Begin: install golang
ENV GOLANG_VERSION 1.23.2
ENV GOLANG_DOWNLOAD_URL https://golang.org/dl/go$GOLANG_VERSION.linux-amd64.tar.gz
ENV GOLANG_SHA256_CHECKSUM 979694c2c25c735755bf26f4f45e19e64e4811d661dd07b8c010f7a8e18adfca
ENV GOPATH $HOME/go
ENV PATH $PATH:/usr/local/go/bin:$GOPATH/bin
RUN curl -fsSL "$GOLANG_DOWNLOAD_URL" -o golang.tar.gz && \
    echo "$GOLANG_SHA256_CHECKSUM golang.tar.gz" | sha256sum -c - && \
    sudo tar -C /usr/local -xzf golang.tar.gz && \
    rm golang.tar.gz && \
    mkdir -p "$GOPATH/src" "$GOPATH/bin" && chmod -R 777 "$GOPATH"
#End: install golang

#Begin: install delve
RUN go install github.com/go-delve/delve/cmd/dlv@v1.23.1
#End: install delve

#Begin: install seqkit
RUN go install github.com/shenwei356/seqkit/v2/seqkit@v2.8.2
#End: install seqkit

WORKDIR $HOME/go/src/github.com/shenwei356/seqkit