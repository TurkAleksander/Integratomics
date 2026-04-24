#!/bin/bash
#You may need to run as sudo if you're experiencing problems with permissions
docker run --rm -v /some/input/dir:/input -v /some/output/dir:/output docker_img_integ /input /output
