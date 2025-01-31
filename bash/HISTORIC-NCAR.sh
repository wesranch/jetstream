#!/bin/bash

set -e
scenario_file="FULLscenario_DGS_hist_ncar.txt"

if [ ! -f "$scenario_file" ]; then
    echo "Error: Scenario file '$scenario_file' does not exist."
    exit 1
fi

#run
dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Console.dll "$scenario_file"
echo "running historic ncar"

# Pause equivalent for user input at the end (optional)
read -p "Press any key to exit..."
