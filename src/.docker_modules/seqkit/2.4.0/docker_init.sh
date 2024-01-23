#!/bin/sh
# docker pull lbmc/seqkit:2.4.0
# docker build src/.docker_modules/seqkit/2.4.0 -t 'lbmc/seqkit:2.4.0'
# docker push lbmc/seqkit:2.4.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/seqkit:2.4.0" --push src/.docker_modules/seqkit/2.4.0

# docker pull xgrand/seqkit:2.4.0
docker build src/.docker_modules/seqkit/2.4.0 -t 'xgrand/seqkit:2.4.0'
docker push xgrand/seqkit:2.4.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/seqkit:2.4.0" --push src/.docker_modules/seqkit/2.4.0