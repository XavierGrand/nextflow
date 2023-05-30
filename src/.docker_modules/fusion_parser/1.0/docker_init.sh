#!/bin/sh
# docker pull xgrand/fusion_parser:1.0
docker build src/.docker_modules/fusion_parser/1.0 -t 'xgrand/fusion_parser:1.0'
docker push xgrand/fusion_parser:1.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/fusion_parser:1.0" --push src/.docker_modules/fusion_parser/1.0