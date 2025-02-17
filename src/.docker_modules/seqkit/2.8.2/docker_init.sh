#!/bin/sh
# docker pull lbmc/seqkit:2.8.2
# docker build src/.docker_modules/seqkit/2.8.2 -t 'lbmc/seqkit:2.8.2'
# docker push lbmc/seqkit:2.8.2
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/seqkit:2.8.2" --push src/.docker_modules/seqkit/2.8.2

# docker pull xgrand/seqkit:2.8.2
docker build src/.docker_modules/seqkit/2.8.2 -t 'xgrand/seqkit:2.8.2'
docker push xgrand/seqkit:2.8.2
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/seqkit:2.8.2" --push src/.docker_modules/seqkit/2.8.2