#!/bin/sh
# docker pull xgrand/viral-track:1.0
docker build src/.docker_modules/viral-track/1.0 -t 'xgrand/viral-track:1.0'
docker push xgrand/viral-track:1.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/viral-track:1.0" --push src/.docker_modules/viral-track/1.0
