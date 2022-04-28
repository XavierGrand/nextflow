#!/bin/sh
docker pull lbmc/htseq:1.99.2
docker build src/.docker_modules/htseq/1.99.2 -t 'lbmc/htseq:1.99.2'
docker push lbmc/htseq:1.99.2
