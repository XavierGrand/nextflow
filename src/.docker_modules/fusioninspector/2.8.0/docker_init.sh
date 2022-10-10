#!/bin/sh
# docker pull xgrand/fusioninspector:2.8.0
docker build src/.docker_modules/fusioninspector/2.8.0 -t 'xgrand/fusioninspector:2.8.0'
docker push xgrand/fusioninspector:2.8.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/fusioninspector:2.8.0" --push src/.docker_modules/fusioninspector/2.8.0
