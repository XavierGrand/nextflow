#!/bin/sh
 docker pull lbmc/minimap2:2.24
 docker build src/.docker_modules/minimap2/2.24 -t 'lbmc/minimap2:2.24'
 docker push lbmc/minimap2:2.24
 docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/minimap2:2.24" --push src/.docker_modules/minimap2/2.24