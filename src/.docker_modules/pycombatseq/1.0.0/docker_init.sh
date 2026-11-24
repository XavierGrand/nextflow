#!/bin/sh
#docker pull xgrand/pycombatseq:1.0.0
docker build src/.docker_modules/pycombatseq/1.0.0 -t 'xgrand/pycombatseq:1.0.0'
docker push xgrand/pycombatseq:1.0.0
#docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/pycombatseq:1.0.0" --push src/.docker_modules/pycombatseq/1.0.0