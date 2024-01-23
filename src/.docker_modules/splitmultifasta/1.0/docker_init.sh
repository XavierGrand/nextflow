#!/bin/sh
# docker pull xgrand/splitmultifasta:1.0
docker build src/.docker_modules/splitmultifasta/1.0 -t 'xgrand/splitmultifasta:1.0'
docker push xgrand/splitmultifasta:1.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/splitmultifasta:1.0" --push src/.docker_modules/splitmultifasta/1.0