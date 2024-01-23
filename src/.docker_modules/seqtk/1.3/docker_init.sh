#!/bin/sh
# docker pull xgrand/seqtk:1.3
docker build src/.docker_modules/seqtk/1.3 -t 'xgrand/seqtk:1.3'
docker push xgrand/seqtk:1.3
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/seqtk:1.3" --push src/.docker_modules/seqtk/1.3