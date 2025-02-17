#!/bin/sh
# docker pull xgrand/cellranger:9.0.0
docker build src/.docker_modules/cellranger/9.0.0 -t 'xgrand/cellranger:9.0.0'
docker push xgrand/cellranger:9.0.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/cellranger:9.0.0" --push src/.docker_modules/cellranger/9.0.0
