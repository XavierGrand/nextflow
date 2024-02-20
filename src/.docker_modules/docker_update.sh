#!/bin/sh
fd "Dockerfile" src/.docke_modules | perl -pe 's|.*docker_modules/(.*)/(.*)/Dockerfile|\1:\2|g' | awk '{system("docker tag "$0" lbmc/" $0)}'
