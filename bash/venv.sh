#!/bin/bash

# working directory and virtual environment path
wd=$(pwd)
venv_path="$wd/env"

# see if python3 is installed
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