#!/bin/sh
# docker pull lbmc/fastp:0.24.0
# docker build src/.docker_modules/fastp/0.24.0 -t 'lbmc/fastp:0.24.0'
# docker push lbmc/fastp:0.24.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.24.0" --push src/.docker_modules/fastp/0.24.0

# docker pull xgrand/fastp:0.24.0
docker build src/.docker_modules/fastp/0.24.0 -t 'xgrand/fastp:0.24.0'
docker push xgrand/fastp:0.24.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.24.0" --push src/.docker_modules/fastp/0.24.0