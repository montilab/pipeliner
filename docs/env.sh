# Create a local copy where
$ mkdir /path/to/conda_envs

# Update this path in environment variables
$ nano ~/.bashrc

# Add the following in file
SCC_CONDA_ENVS=/path/to/conda_envs
export SCC_CONDA_ENVS

# After saving
$ source ~/.bashrc

# Comfirm the path is correct
$ echo $SCC_CONDA_ENVS

# Go to your new local conda envs when creating or activating new environments
$ cd $SCC_CONDA_ENVS

# Create new local environment call pipeliner
$ conda create --prefix=$SCC_CONDA_ENVS/pipeliner python=2.7

# Activate environment
$ source activate pipeliner

# Install pipeliner dependencies
$ conda install --prefix=$SCC_CONDA_ENVS/pipeliner --channel "anfederico" trim-galore fastqc star multiqc samtools rseqc stringtie hisat2 htseq subread numpy pandas bioconductor-biobase

# Check to make sure dependencies are installed correctly
$ conda list --prefix=$SCC_CONDA_ENVS
