#!/bin/sh
# docker pull xgrand/samtools:1.20
docker build src/.docker_modules/samtools/1.20 -t 'xgrand/samtools:1.20'
docker push xgrand/samtools:1.20
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/samtools:1.30" --push src/.docker_modules/samtools/1.20
