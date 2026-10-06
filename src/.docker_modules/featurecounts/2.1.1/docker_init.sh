#!/bin/sh
# docker pull xgrand/featurecounts:2.1.1
docker build src/.docker_modules/featurecounts/2.1.1 -t 'xgrand/featurecounts:2.1.1'
docker push xgrand/featurecounts:2.1.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/featurecounts:2.1.1" --push src/.docker_modules/featurecounts/2.1.1
