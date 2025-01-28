#!/bin/bash

# run the container
sudo docker container run -d -p 127.0.0.1:8080:7681 -v ~/landis/input:/container/data esiil/landis_v8:latest

# output container name
container_name=$(sudo docker container ls --format '{{.Names}}' | tail -n 1)
echo "Container started: $container_name"

#user defines scenarios and reps
#read -p "Which scenarios do you want to run (full basename Comma Separated): " batch_files
read -p "Number of replicates to run?" num_replicates
read -p "access ID?:" access_id

export batch_files
export access_id
export num_replicates

echo "scenarios: "$batch_files""
echo "reps: "$num_replicates""

#install python boto3 for the python script
pip3 install boto3 #not working

mv ~run-landis.py ~landis/input/run-landis.py
python3 ~landis/input/run-landis.py \
    --container_name "$container_name" \
    --batch_files "$batch_files" \
    --access_id "$access_id" \
    --num_replicates "$num_replicates"
