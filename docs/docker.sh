# Find list of docker repositories
docker images

# Using the miniconda fresh docker image
docker pull continuumio/miniconda
docker run -i -t continuumio/miniconda /bin/bash

# But you could also just do the following for any repository
docker run -i -t <repository> /bin/bash

# Docker image is now activating, can access conda environment from here
conda list

# Install your packages
conda install --channel "anfederico" trim-galore fastqc star multiqc samtools rseqc stringtie hisat2 htseq subread numpy pandas bioconductor-biobase

# while in the docker image, grab the id (looks like c16378f943fe)
# Create a new repository with updated conda environment (copy)
docker commit c16378f943fe <repository>

docker login
docker images
docker tag bb38976d03cf anfederico/pipeliner:<repository>
docker push anfederico/pipeliner

# To download on new computer
docker pull anfederico/<repository>
