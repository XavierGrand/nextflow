#!/bin/sh
docker pull xgrand/htseq:2.0.5
docker build src/.docker_modules/htseq/2.0.5 -t 'xgrand/htseq:2.0.5'
docker push xgrand/htseq:2.0.5
