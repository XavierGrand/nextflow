#!/bin/sh
docker pull lbmc/guppy-cpu:5.0.11
docker build src/.docker_modules/guppy-cpu/5.0.11 -t 'lbmc/guppy-cpu:5.0.11'
docker push lbmc/guppy-cpu:5.0.11
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/guppy-cpu:5.0.11" --push src/.docker_modules/guppy-cpu/5.0.11

