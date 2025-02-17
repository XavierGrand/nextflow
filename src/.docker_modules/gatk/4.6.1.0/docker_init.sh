#!/bin/sh
# docker pull xgrand/gatk:4.6.1.0
docker build src/.docker_modules/gatk/4.6.1.0 -t 'xgrand/gatk:4.6.1.0'
docker push xgrand/gatk:4.6.1.0
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/gatk:4.6.1.0" --push src/.docker_modules/gatk/4.6.1.0
