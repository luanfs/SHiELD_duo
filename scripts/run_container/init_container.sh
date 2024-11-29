#!/bin/bash
# Bash script to initialize the container

# get directory
cd ..
container_dir=$(pwd)
cd -

# Initialize container
docker run --privileged --cap-add=SYS_PTRACE --ulimit='stack=-1:-1' -v $container_dir/SHiELD_SRC:/SHiELD_SRC -it gfdlfv3/shield-dev bash
