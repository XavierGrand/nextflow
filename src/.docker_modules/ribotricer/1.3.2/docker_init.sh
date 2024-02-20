#!/bin/sh
docker pull lbmc/ribotricer:1.3.2
# docker build src/.docker_modules/ribotricer/1.3.2 -t 'lbmc/ribotricer:1.3.2'
# docker push lbmc/ribotricer:1.3.2
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/ribotricer:1.3.2" --push src/.docker_modules/ribotricer/1.3.2
