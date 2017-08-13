# Various scripts for Ubuntu 16.04 environment setup

# First install miniconda2
https://conda.io/miniconda.html

# Edit the path so miniconda3 is invoked when calling 'conda'
gedit ~/.bashrc

# Create a new environment
cd ~/miniconda2/envs
conda create -n pipeliner python=2.7

# Make sure you see it in the correct location
conda info --envs

# Activate it
source activate pipeliner

# Install dependencies
conda install -c bioconda java-jdk fastqc trim-galore star multiqc samtools rseqc stringtie

# Make sure Java 7/8 is installed
java -version

# Install Nextflow executable
cd ~/Documents/rnaseq/pipeliner/Gallus_example
curl -s https://get.nextflow.io | bash

# Test Nextflow executable
./nextflow run hello

