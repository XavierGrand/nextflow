#!/bin/sh
docker pull lbmc/minimap2:2.20
docker build src/.docker_modules/minimap2/2.20 -t 'lbmc/minimap2:2.20'
docker push lbmc/minimap2:2.20
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/minimap2:2.20" --push src/.docker_modules/minimap2/2.20