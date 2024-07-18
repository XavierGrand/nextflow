#!/bin/sh
# docker pull xgrand/ksnp:4.1
docker build src/.docker_modules/ksnp/4.1 -t 'xgrand/ksnp:4.1'
docker push xgrand/ksnp:4.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/ksnp:4.1" --push src/.docker_modules/ksnp/4.1