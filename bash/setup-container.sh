#!/bin/bash

# run the container
sudo docker container run -d -p 127.0.0.1:8080:7681 -v ~/landis/input:/home/user esiil/landis_v8:latest

# output container name
container_name=$(sudo docker container ls --format '{{.Names}}' | tail -n 1)
echo "Container started: $container_name"
#access username
read -p "access ID/username?:" access_id
export access_id

#install python boto3 for the python script
pip3 install boto3
python3 ~landis/input/run-landis.py \
    --container_name "$container_name" \
	--access_id "$access_id"