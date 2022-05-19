#!/bin/sh
# docker pull lbmc/stringtie:2.1.1
docker build src/.docker_modules/stringtie/2.1.0 -t 'lbmc/stringtie:2.1.1'
docker push lbmc/stringtie:2.1.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/stringtie:2.1.1" --push src/.docker_modules/stringtie/2.1.1