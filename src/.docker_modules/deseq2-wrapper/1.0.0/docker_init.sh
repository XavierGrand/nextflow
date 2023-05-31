#!/bin/sh
# docker pull xgrand/deseq2-wrapper:1.0.0
docker build src/.docker_modules/deseq2-wrapper/1.0.0 -t 'xgrand/deseq2-wrapper:1.0.0'
docker push xgrand/deseq2-wrapper:1.0.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/deseq2-wrapper:1.0.0" --push src/.docker_modules/deseq2-wrapper/1.0.0