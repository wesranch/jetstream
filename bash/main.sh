#!/bin/bash

cd /home/wancher/Documents/thesis/bash
./connect-jetstream.sh

sleep 5

./setup-jetstream
./setup-container

#change port on setup container and run again