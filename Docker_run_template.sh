#!/bin/bash
#You may need to run as sudo if you're experiencing problems with permissions
docker run --rm -v /scratch/PROJECTS/KT_obesity_integ/For_integ:/input -v /scratch/PROJECTS/KT_obesity_integ/Integ_output:/output docker_img_integ /input /output