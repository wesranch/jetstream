#!/bin/bash
set -euo pipefail

trap 'echo "Error occurred at line $LINENO"; exit 1;' ERR

# Define working directory and virtual environment path
wd=$(pwd)
venv_path="$wd/env"

# Ensure python3 is installed
if ! command -v python3 &> /dev/null; then
  echo "python3 is required but not installed. Please install it."
  exit 1
fi

# check if the virtual environment already exists
if [ -d "$venv_path" ]; then
  echo "Activating virtual environment..."
  source "$venv_path/bin/activate"
else
  echo "Creating virtual environment at: $venv_path"
  python3 -m venv env
  source "$venv_path/bin/activate"
fi

#environment variables
export AWS_ACCESS_KEY_ID="your-access-key-id"
export AWS_SECRET_ACCESS_KEY="your-secret-access-key"
export S3_ENDPOINT_URL=${S3_ENDPOINT_URL:-"https://usgs2.osn.mghpcc.org"}
export SSH_KEYS_DIR="$HOME/.ssh"

#params for aws
# read -p "s3 bucket name: " bucket_name
# read -p "where do you want your training data inside landis-ii-input bucket?: " training_data_folder
# read -p "where do you want your random forest CODE inside landis-ii-input bucket?: " bucket_code_folder
# read -p "where are your rasters needed for predictions inside landis-ii-input bucket?: " bucket_raster_folder

# #user defines paths to files
# read -p "relative path to your training data on your local machine?: " training_data_path
# echo "uploading most recent training data to s3..."

# aws s3 sync "$training_data_path" s3://"$bucket_name"/"$training_data_folder" --endpoint-url "$S3_ENDPOINT_URL"

# read -p "Provide the relative path to your random forest code: " scripts_path
# echo "uploading most recent model script to s3..."
# aws s3 sync "$scripts_path" s3://"$bucket_name"/"$bucket_code_folder" --endpoint-url "$S3_ENDPOINT_URL"

#SSH into jetstream
read -p "Enter your ACCESS username: " access_id
read -p "Paste the public IP address of your jetstream deployment: " pub_ip

export access_id
export bucket_name
#export bucket_code_folder
#export training_data_folder
#export bucket_raster_folder

ssh -i "$SSH_KEYS_DIR" -o SendEnv=AWS_ACCESS_KEY_ID -o SendEnv=AWS_SECRET_ACCESS_KEY \
  -o SendEnv=bucket_name -o SendEnv=S3_ENDPOINT_URL "$access_id"@"$pub_ip" << 'ENDSSH'

  sleep 10
  
  #see if AWS is already installed
  if ! command -v aws &> /dev/null
  then
    echo "AWS CLI not found. Installing..."
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    unzip awscliv2.zip
    sudo ./aws/install
  else
    echo "AWS CLI already installed."
  fi

  # Check if AWS credentials are already configured
  #unexpected EOF
  if [ ! -f "$HOME/.aws/credentials" ]; then
    echo "AWS credentials not configured. Configuring now..."
    mkdir -p ~/.aws
    cat <<'EOF' > ~/.aws/credentials
    [default]
    aws_access_key_id=$AWS_ACCESS_KEY_ID
    aws_secret_access_key=$AWS_SECRET_ACCESS_KEY
    region=us-west-2
    output=json
    EOF
    echo "AWS credentials stored in ~/.aws/credentials"
  else
    echo "AWS credentials already configured."
  fi

  echo "downloading training data from s3..."
  mkdir -p training-data
  aws s3 cp s3://"$bucket_name"/rs-data/training-data/ --endpoint-url "$S3_ENDPOINT_URL" training-data --recursive

  echo "downloading rasters to predict on from s3..."
  mkdir -p rasters
  aws s3 cp s3://"$bucket_name"/rs-data/rasters/ --endpoint-url "$S3_ENDPOINT_URL" rasters --recursive

  #Location of R scripts
  mkdir -p scripts
  aws s3 cp s3://"$bucket_name"/scripts/ --endpoint-url "$S3_ENDPOINT_URL" scripts --recursive

  echo "creating model output directory at models/"
  mkdir -p models

  echo "creating folder for output graphs at figures/"
  mkdir -p figures

  # echo "Running random forest classification code..."
  # Rscript scripts/random-forest.R
  # echo "Random forest code completed at: timestamp"

  # echo "models saved to "models/""
  # echo "model performance figures saved to "figures/""


  # echo "Running random forest prediction code..."
  # Rscript scripts/random-forest-predict.R
  # echo "prediction maps saved to "rasters/""

ENDSSH

#make it executable
#chmod +x rf.sh
