#!/bin/bash

#install wget
sudo apt update
sudo apt install -y curl unzip wget

#aws
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
aws --version
aws configure

#gdrive
sudo apt install brew
brew install gdrive
gdrive account add #client id on google cloud profile

#got these numbers from the share url
gdrive files list --query "'1TUtmF2tGdFFjvt2xMB2j-gfhHnTTRvVW' in parents"

#download the files on the VM
mkdir -p rs-data/tiles rs-data/mosaics
for FILE_ID in $(gdrive files list --query "'1TUtmF2tGdFFjvt2xMB2j-gfhHnTTRvVW' in parents" --skip-header | awk '{print $1}'); do
    FILE_NAME=$(gdrive files info $FILE_ID | grep "Name" | awk -F': ' '{print $2}')
    wget --no-check-certificate "https://drive.google.com/uc?export=download&id=$FILE_ID" -O "/home/wrancher/rs-data/tiles/$FILE_NAME"
done

#use this if not using Rstudio image
#module avail R
#module load R/4.2.1
#Rscript mosaic_rasters.R