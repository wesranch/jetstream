#!/bin/bash
read -p "image name:" image_name
docker build --platform linux/amd64 -t "$image_name" .  
export IMAGE_NAME="$image_name"
