#!/bin/sh
docker pull xgrand/star:2.7.8a
# docker build src/.docker_modules/star/2.7.8a/ -t 'xgrand/star:2.7.8a'
# docker push xgrand/star:2.7.8a
docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/star:2.7.8a" --push src/.docker_modules/star/2.7.8a