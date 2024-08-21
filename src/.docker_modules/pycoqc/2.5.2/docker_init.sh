#!/bin/sh
#docker pull xgrand/pycoqc:2.5.2
docker build src/.docker_modules/pycoqc/2.5.2 -t 'xgrand/pycoqc:2.5.2'
docker push xgrand/pycoqc:2.5.2
#docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/pycoqc:2.5.2" --push src/.docker_modules/pycoqc/2.5.2