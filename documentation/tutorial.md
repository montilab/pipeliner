## Nextflow Tutorial
#### Pre-requisites
Pipeliner requires Java, Nextflow, and Anaconda for implementation. All other tools for implementation are wrapped in the conda environment described below. 

*1. Test Nextflow*

Make sure you have Java 7/8 installed and then install Nextflow to any working directory. Test the Nextflow executable before continuuing. This is just to make sure your environment is compatible with a nextflow executable. You will download another one later in the tutorial.
```bash
java -version
cd path/to/wd
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

*2. Download Conda*

- Local Machine  
Conda is available through [Anaconda](https://www.continuum.io/downloads). Download the distribution pre-packaged with Python 2.7. If this is your first time working with conda, you may need to edit your configuration paths to ensure anaconda is invoked when calling `conda`

- Shared Computing Cluster
Enable `conda` by loading a pre-installed version of Anaconda with `module load anaconda2`. This will load the latest pre-installed version pre-packaged with Python 2.7.

#### Running the pipeline

*3. Clone Tutorial Repository*

```bash
git clone https://github.com/montilab/pipeliner
cd pipeliner/RNA-seq
```

*4. Create Conda Environment*

```bash
conda env create -f envs/linux_env.yml 
source activate pipeliner

or

conda env create -f envs/osx_env.yml 
source activate pipeliner
```

*5. Run Workflow*
As of now, you must explicitly declare a few input paths in the following files:  
*pipeliner/RNA-seq/nextflow.config*  
*pipeliner/RNA-seq/ggal_data/ggal_reads.csv*

Once finished, run the pipeline

```bash
curl -s https://get.nextflow.io | bash
./nextflow main.nf -c nextflow.config
```
