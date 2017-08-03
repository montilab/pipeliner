# Various scripts for Ubuntu 16.04 environment setup

# First install miniconda3
https://conda.io/miniconda.html

# Edit the path so miniconda3 is invoked when calling 'conda'
gedit ~/.bashrc

# Create a new environment
cd ~/miniconda3/envs
conda create -n pypliner3 python=3.5

# Make sure you see it in the correct location
conda info --envs

# Activate it
source activate pypliner3

# Install dependencies
conda install -c bioconda java-jdk
conda install -c bioconda fastqc
conda install -c bioconda trim-galore
conda install -c bioconda star
conda install -c bioconda multiqc

# Make sure Java 7/8 is installed
java -version

# Install Nextflow executable
cd ~/Documents/pythonvillage/pypliner3/Gallus_example
curl -s https://get.nextflow.io | bash

# Test Nextflow executable
./nextflow run hello

