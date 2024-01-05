#!/bin/sh
# docker pull lbmc/readthroughbincounts:0.1
docker build src/.docker_modules/readthroughbincounts/0.1 -t 'lbmc/readthroughbincounts:0.1'
# docker push lbmc/readthroughbincounts:0.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "lbmc/readthroughbincounts:0.1" --push src/.docker_modules/readthroughbincounts/0.1

# docker pull xgrand/readthroughbincounts:0.1
# docker build src/.docker_modules/readthroughbincounts/0.1 -t 'xgrand/readthroughbincounts:0.1'
# docker push xgrand/readthroughbincounts:0.1
# docker buildx build --platform linux/amd64,linux/arm64 -t "xgrand/readthroughbincounts:0.1" --push src/.docker_modules/readthroughbincounts/0.1
