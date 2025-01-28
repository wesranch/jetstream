#!/bin/bash

# Make sure the script exits on error
set -e

# Add the alias for landis-ii-8 if not already added
if ! grep -q "alias landis-ii-8" ~/.bashrc; then
    echo "alias landis-ii-8='../../../opt/Core-Model-v8-LINUX-main/build/Release/Landis.Console'" >> ~/.bashrc
    # Reload bashrc
    source ~/.bashrc
fi

# Run Landis-II with the specified scenario file
scenario_file="FULLscenario_DGS_fut_gfdl.txt"

if [ ! -f "$scenario_file" ]; then
    echo "Error: Scenario file '$scenario_file' does not exist."
    exit 1
fi

# Call the Landis-II alias
landis-ii-8 "$scenario_file"
echo "running future gfdl"
# Pause equivalent for user input at the end (optional)
read -p "Press any key to exit..."
