#!/bin/sh
docker pull mlepetit/saanity
docker build src/.docker_modules/mlepetit/sanity -t 'lbmc/sanity'
docker push lbmc/sanity
