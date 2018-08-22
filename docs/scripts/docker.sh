# Find list of docker repositories
docker images

# Using the miniconda fresh docker image
docker pull continuumio/miniconda
docker run -i -t continuumio/miniconda /bin/bash

# Docker image is now activating, can access conda environment from here
conda list

# Install your packages
conda install --channel "anfederico" trim-galore fastqc star multiqc samtools rseqc

# Clone an entire environment
cd root
git clone https://github.com/montilab/pipeliner
/root/pipeliner/envs/osx_env.yml
cd /
conda env update -f=/root/pipeliner/envs/linux_env.yml

# Setup pipeliner paths and download Nextflow executable
source activate pipeliner
python pipeliner/scripts/paths.py
cd pipeliner/pipelines
curl -s https://get.nextflow.io | bash

# while in the docker image, grab the id (looks like c16378f943fe)
# Create a new repository with updated conda environment (copy)
docker commit c16378f943fe latest

# Push to docker hub
docker login
docker images
docker tag bb38976d03cf pipeliner/images:latest
docker push pipeliner/images

# To download on new computer
docker pull pipeliner/images
docker run -i -t pipeliner/images /bin/bash

# Pipeliner with all dependencies installed
source activate pipeliner
cd /root/pipeliner
