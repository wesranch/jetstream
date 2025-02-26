#!/bin/bash

set -e

mkdir -p mosaics
mkdir -p tiles

# use sudo if writing for container
if ! command -v brew &> /dev/null; then
    echo "Installing Brew"
    sudo apt update
    sudo apt-get install build-essential
    sudo apt install git -y
    sudo apt-get install tar -y
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
    brew --version
else
    echo "brew already installed..."
fi

# install Google Drive and Gdownload libs
if ! command -v gdrive &> /dev/null; then
    brew install gdrive
else 
    echo "gdrive already installed..."
fi

if ! command -v pipx &> /dev/null; then
    brew install pipx
else
    echo "pipx already installed..."
fi

if ! pipx run gdown --version &> /dev/null; then
    pipx install gdown
else
    echo "gdown already installed..."
fi

# set up Google Drive account
# for running on remote server see README: https://github.com/glotlabs/gdrive
if ! gdrive about &> /dev/null; then
    echo "No account found. Adding account..."
    echo "Client ID and secret are found on Google Cloud/ API and Services/ Credentials (Click Auth for Secret Key)"
    
    #export GOOGLE_APP_CREDS="/home/$access_id/gdrive_export-wesley_rancher_2023_owu_edu.tar"
    #tar -xf "$GOOGLE_APP_CREDS" -C /home/$access_id/
    gdrive account add --service-account "$GOOGLE_APP_CREDS"
    #gdrive account add
else
    echo "Google Drive account added..."
fi

# list files in account and download them // Works fastest with archived files
read -p "Desired destination to download to (full path): " DEST_DIR
gdrive files list
read -p "Do you need specific files from a folder? " response

export DEST_DIR
export response

if [[ "$response" == "y" || "$response" == "Y" ]]; then
    read -p "Paste parent folder Id: " FOLDER_ID
    gdrive files list --parent="$FOLDER_ID"


    read -p "Do you want to download all files matching a pattern? (y/n): " pattern_response
    if [[ "$pattern_response" == "y" || "$pattern_response" == "Y" ]]; then
        read -p "Enter the pattern to match: " FILE_PATTERN
        MATCHED_FILES=$(gdrive files list --parent="$FOLDER_ID" | grep "$FILE_PATTERN" | awk '{print $1}')

        if [[ -z "$MATCHED_FILES" ]]; then
            echo "No files found matching the pattern '$FILE_PATTERN'."
            exit 1
        fi

        for FILE_ID in $MATCHED_FILES; do
            gdrive files download "$FILE_ID" --destination="$DEST_DIR"
        done
    else
        read -p "Paste file Ids from parent folder (comma-separated): " FILE_IDS_NEEDED
        IFS=',' read -r -a FILE_IDS_ARRAY <<< "$FILE_IDS_NEEDED"

        for FILE_ID in "${FILE_IDS_ARRAY[@]}"; do
            gdrive files download "$FILE_ID" --destination="$DEST_DIR"
        done
    fi

elif [[ "$response" == "n" || "$response" == "N" ]]; then
    read -p "Paste folder or file Id needed: " FOLDER_OR_FILE
    FILE_INFO=$(gdrive files info "$FOLDER_OR_FILE")

    if echo "$FILE_INFO" | grep -q "Mime: application/vnd.google-apps.folder"; then
        read -p "Download folder recursively? (y/n): " recursive
        if [[ "$recursive" == "y" || "$recursive" == "Y" ]]; then
            gdrive files download --recursive "$FOLDER_OR_FILE" --destination="$DEST_DIR"
        else
            gdrive files download "$FOLDER_OR_FILE" --destination="$DEST_DIR"
        fi
    else
        gdrive files download "$FOLDER_OR_FILE" --destination="$DEST_DIR"
    fi
fi

# arguments to parse into R script
read -p "Unique pattern consistent across all files: " PATTERN
read -p "Output directory for mosaiced images (full-path): " OUT_DIR
read -p "Start year: " START_YEAR
read -p "End year: " END_YEAR

export PATTERN
export OUT_DIR
export START_YEAR
export END_YEAR

# R and dependencies
sudo apt install r-base
R --version
Rscript -e "install.packages("terra")"

echo "Mosaicing tiles..."

# call the mosaic rasters script
Rscript /home/wrancher/mosaic-rasters.R \
     --dest_dir "$DEST_DIR" \
     --pattern "$PATTERN" \
     --start_year "$START_YEAR" \
     --end_year "$END_YEAR" \
     --out_dir "$OUT_DIR"

echo "Mosaicing completed..."

#clean sampling csvs
#Rscript /Users/wancher/Documents/thesis/scripts/cleanSampling.R \
#    --dest_dir "$DEST_DIR" \
#    --pattern "$PATTERN" \
#    --out_dir "$OUT_DIR"

