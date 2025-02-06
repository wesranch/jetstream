#!/bin/bash

#replace with abs path to input files or just do docker run -t IMAGE_NAME
docker run --platform linux/amd64 -v /path/to/input/files:/input -t IMAGE_NAME
