#!/bin/sh
# docker pull lbmc/fastp:0.23.2
# docker build src/.docker_modules/fastp/0.23.2 -t 'lbmc/fastp:0.23.2'
# docker push lbmc/fastp:0.23.2
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.23.2" --push src/.docker_modules/fastp/0.23.2

# docker pull xgrand/fastp:0.23.2
docker build src/.docker_modules/fastp/0.23.2 -t 'xgrand/fastp:0.23.2'
docker push xgrand/fastp:0.23.2
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/fastp:0.23.2" --push src/.docker_modules/fastp/0.23.2