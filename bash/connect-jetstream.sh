#!/bin/bash

# grab local env vars for aws (not working)
#export AWS_ACCESS_KEY_ID="your-access-key-id"
#export AWS_SECRET_ACCESS_KEY="your-secret-access-key"
#export S3_ENDPOINT_URL=${S3_ENDPOINT_URL:-"https://usgs2.osn.mghpcc.org"}

# private ssh in (you may need to cp your .pub key onto JS2 web terminal)
export SSH_KEYS_DIR="$HOME/.ssh/id_rsa"

read -p "Enter ACCESS username: " access_id
read -p "Enter Jetstream VM Public IP: " pub_ip

# secure copy scripts to vm
scp -i "$SSH_KEYS_DIR" setup-jetstream.sh "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" setup-container.sh "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" climate.sh "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" mosaic-gee.sh "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" s3-upload.sh "$access_id@$pub_ip:~/"

#R and python scripts
scp -i "$SSH_KEYS_DIR" ../scripts/climate-data-processing.Rmd "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" ../scripts/random-forest.R "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" ../scripts/run-landis.py "$access_id@$pub_ip:~/"
scp -i "$SSH_KEYS_DIR" ../scripts/mosaic_rasters.R "$access_id@$pub_ip:~/"


#connect and make scripts executable
ssh -i "$SSH_KEYS_DIR" "$access_id@$pub_ip"
chmod +x setup-jetstream.sh
chmod +x setup-container.sh
chmod +x climate.sh
chmod +x mosaic-gee.sh
chmod +x s3-upload.sh