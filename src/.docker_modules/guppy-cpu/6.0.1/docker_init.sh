#!/bin/sh
docker pull lbmc/guppy-cpu:6.0.1
docker build src/.docker_modules/guppy-cpu/6.0.1 -t 'lbmc/guppy-cpu:6.0.1'
docker push lbmc/guppy-cpu:6.0.1
docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/guppy-cpu:6.0.1" --push src/.docker_modules/guppy-cpu/6.0.1
