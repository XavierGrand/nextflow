#!/bin/sh
docker pull xgrand/rtranger:1.0
# docker build src/.docker_modules/rtranger/1.0 -t 'xgrand/rtranger:1.0'
# docker push xgrand/rtranger:1.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/rtrange:1.0" --push src/.docker_modules/rtranger/1.0