#!/bin/bash
docker build -f docker/Dockerfile -t bonesim:latest .
docker run --rm \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/pictures:/workspace/pictures \
    -v $(pwd)/report:/workspace/report \
    -v $(pwd)/calculix:/workspace/calculix \
    bonesim:latest "$@"