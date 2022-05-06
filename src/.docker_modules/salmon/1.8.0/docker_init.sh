#!/bin/sh
docker pull lbmc/salmon:1.8.0
docker build src/.docker_modules/salmon/1.8.0 -t 'lbmc/salmon:1.8.0'
docker push lbmc/salmon:1.8.0
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/salmon:1.8.0" --push src/.docker_modules/salmon/1.8.0