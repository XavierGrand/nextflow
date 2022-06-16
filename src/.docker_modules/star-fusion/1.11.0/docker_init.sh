#!/bin/sh
# docker pull xgrand/star-fusion:1.11.0
docker build src/.docker_modules/star-fusion/1.11.0 -t 'xgrand/star-fusion:1.11.0'
docker push xgrand/star-fusion:1.11.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/star-fusion:1.11.0" --push src/.docker_modules/star-fusion/1.11.0
