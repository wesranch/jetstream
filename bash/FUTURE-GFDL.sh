#!/bin/bash

set -e
scenario_file="FULLscenario_DGS_fut_gfdl.txt"

if [ ! -f "$scenario_file" ]; then
    echo "Error: Scenario file '$scenario_file' does not exist."
    exit 1
fi

#run
dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Console.dll "$scenario_file"
echo "running future gfdl"

# Pause equivalent for user input at the end (optional)
read -p "Press any key to exit..."
