#!/bin/sh
docker pull mlepetit/sanity
docker build src/.docker_modules/sanity -t 'lbmc/sanity'
docker push lbmc/sanity
