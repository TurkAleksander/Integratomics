#!/bin/bash
#You may need to run as sudo if you're experiencing problems with permissions
docker run --rm -v /path/to/local/input:/input -v /path/to/local/output:/output docker_img_integ /input /output