#!/bin/bash
#
# This is the docker build command. Feel free to modify any tags here
docker build --no-cache --network=host -t spotlight_dnv:spotlight .
# Note: You might have to run this with admin privilages (sudo)
