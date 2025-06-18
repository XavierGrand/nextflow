#!/bin/sh
docker pull xgrand/htseq:2.0.9
docker build src/.docker_modules/htseq/2.0.9 -t 'xgrand/htseq:2.0.9'
docker push xgrand/htseq:2.0.9
