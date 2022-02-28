#!/bin/sh
docker pull lbmc/hisat2:2.2.1
docker build src/.docker_modules/hisat2/2.2.1 -t 'lbmc/hisat2:2.2.1'
docker push lbmc/hisat2:2.2.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/hisat2:2.2.1" --push src/.docker_modules/hisat2/2.2.1
