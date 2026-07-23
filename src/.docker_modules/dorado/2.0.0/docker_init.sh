#!/bin/sh
# docker pull xgrand/dorado:2.0.0
docker build src/.docker_modules/dorado/2.0.0 -t 'xgrand/dorado:2.0.0'
docker push xgrand/dorado:2.0.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/dorado:2.0.0" --push src/.docker_modules/dorado/2.0.0
