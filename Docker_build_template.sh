#!/bin/bash
#Specify the directory where the Dockerfile is located
cd /path/to/dockerfile_directory
#Build the Docker image, don't forget to specify the name of the image here and later when running the container
#The dot (.) at the end of the command is important, it tells Docker to look for the Dockerfile in the current directory
#Added build argument so the Dockerfile always downloads newest Git repository
docker build --build-arg CACHEBUST=$(date +%s) -t docker_img_integ -f Dockerfile_integration .
