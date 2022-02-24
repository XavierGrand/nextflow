#!/bin/sh
docker pull lbmc/bedtools:2.30.0
docker build src/.docker_modules/bedtools/2.30.0 -t 'lbmc/bedtools:2.30.0'
docker push lbmc/bedtools:2.30.0
