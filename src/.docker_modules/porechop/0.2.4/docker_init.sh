#!/bin/sh
docker pull lbmc/porechop:0.2.4
docker build src/.docker_modules/porechop/0.2.4 -t 'lbmc/porechop:0.2.4'
docker push lbmc/porechop:0.2.4
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/porechop:0.2.4" --push src/.docker_modules/porechop/0.2.4