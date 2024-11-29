#!/bin/bash
# Bash script to runt the sw container

# get directory
cd ..
container_dir=$(pwd)
cd -

# Initialize container
docker run --privileged --cap-add=SYS_PTRACE --ulimit stack=51200000:102400000 -v $container_dir/SHiELD_SRC:/SHiELD_SRC -it sw bash

