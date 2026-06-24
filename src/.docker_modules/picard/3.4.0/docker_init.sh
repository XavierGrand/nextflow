#!/bin/sh
# docker pull xgrand/picard:3.4.0
docker build src/.docker_modules/picard/3.4.0 -t 'xgrand/picard:3.4.0'
docker push xgrand/picard:3.4.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/picard:3.4.0" --push src/.docker_modules/picard/2.18.11
