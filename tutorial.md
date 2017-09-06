## Running the pipeline

#### Pre-requisites

Pipeliner requires Java, Nextflow, and Anaconda for implementation. All other tools for implementation are wrapped in the environment described below. 

*1. Download Nextflow*

Make sure you have Java 7/8 and install Nextflow to a working directory. Test the Nextflow executable before continuuing.
```bash
java -version
cd path/to/wd
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

*2. Download Conda*

`Conda` is available through `Anaconda <https://www.continuum.io/downloads>` and  `Miniconda <https://conda.io/miniconda.html>`

If you're using a module system such as on the shared computing cluster (SCC) at Boston University you can just load a preinstalled version::

   module purge
   module load anaconda2/4.3.0


*2. Create Conda Environment (rna_env)*

To install a basic development environment, download rna_env from Pipeliner repository in anaconda cloud::

  conda env create pipeliner/rna_env
  source activate rna_env
  
There is also an option to create a './rna_env folder and store all the files needed to run the pipeline as such::

  conda create \
    -p ./rna_env \
    -c https://anaconda.org/Pipeliner \
    --yes \


Running the pipeline
====================

Activate the environment (follow instructions above to create environment)::
 
  source activate rna_env

OR::

  source activate ./rna_env
  
Pipeliner consists of a main nextflow script parametrized using a configuration file. The configuration file includes all parameters necessary to run the pipeline including  parameters to direct the path of files and results, as well as selecting specific tools and processes to run. The example *Gallus gallus* (chicken) dataset is all within the ggal folder. All files necessary to run this example is in the folder "Gallus_Example".


Once the appropriate paths and tools have been set in the configuration file, you can run the pipeline with the configuration file::

  nextflow main.nf -c config

If an issue arises and pipline stops running, you can check the nextflow log file that is created during each run to determine the problem. Note: A new log file is created after each run and the most recent log file is not numbered::

  vim .nextflow.log #recent log file

Once you have resolved the issue, the pipeline can be resumed using the same command to run the pipeline with '-resume' argument::

  nexflow main.nf -c ./Gallus_Example/config_ggal -resume 

For more information, check out the manual provided.. 
