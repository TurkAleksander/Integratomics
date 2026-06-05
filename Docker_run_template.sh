#!/bin/bash
#You may need to run as sudo if you're experiencing problems with permissions
docker run --rm -v /some/input/dir:/input -v /some/weights/dir:/weights -v /some/output/dir:/output docker_img_integ_weighted /input /weights /output
