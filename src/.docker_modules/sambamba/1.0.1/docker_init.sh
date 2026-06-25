#!/bin/sh
# docker pull xgrand/sambamba:1.0.1
docker build src/.docker_modules/sambamba/1.0.1 -t 'xgrand/sambamba:1.0.1'
docker push xgrand/sambamba:1.0.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/sambamba:1.0.1" --push src/.docker_modules/sambamba/1.0.1
