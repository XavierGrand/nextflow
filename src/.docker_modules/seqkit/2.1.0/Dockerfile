FROM ubuntu:16.04

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
ENV GOLANG_DOWNLOAD_URL https://golang.org/dl/go1.17.9.linux-amd64.tar.gz
ENV GOPATH $HOME/go
ENV PATH $PATH:/usr/local/go/bin:$GOPATH/bin
RUN curl -fsSL "$GOLANG_DOWNLOAD_URL" -o golang.tar.gz && \
    sudo tar -C /usr/local -xzf golang.tar.gz && \
    rm golang.tar.gz && \
    mkdir -p "$GOPATH/src" "$GOPATH/bin" && chmod -R 777 "$GOPATH"
#End: install golang

#Begin: install delve
RUN go install github.com/go-delve/delve/cmd/dlv@latest 
#1.8.3
#End: install delve

#Begin: install seqkit
RUN go install github.com/shenwei356/seqkit/v2/seqkit@v2.1.0

#2.1.0
#End: install seqkit

WORKDIR $HOME/go/src/github.com/shenwei356/seqkit