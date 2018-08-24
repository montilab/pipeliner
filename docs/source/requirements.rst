.. _requirements-page:

Requirements
============

.. toctree::
   :maxdepth: 3
   :caption: Contents:

The `Pipeliner` framework requires `Nextflow` and `Anaconda`. `Nextflow` can be used on any POSIX compatible system (Linux, OS X, etc). It requires BASH and `Java 8 (or higher) <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_ to be installed. Third-party software tools used by individual pipelines will be installed and managed through a Conda virtual environment.

Testing Nextflow
----------------
Before continuuing, test to make sure your environment is compatible with a Nextflow executable. 

.. note:: You will download another one later when you clone the repository

Make sure your Java installation is version 8 or higher::

   java -version

Create a new directory and install/test `Nextflow`::

   mkdir nf-test
   cd nf-test
   curl -s https://get.nextflow.io | bash
   ./nextflow run hello

Output::

   N E X T F L O W  ~  version 0.31.0
   Launching `nextflow-io/hello` [sad_curran] - revision: d4c9ea84de [master]
   [warm up] executor > local
   [4d/479eec] Submitted process > sayHello (4)
   [a8/4bc038] Submitted process > sayHello (2)
   [17/5be64e] Submitted process > sayHello (3)
   [ee/0d879f] Submitted process > sayHello (1)
   Hola world!
   Ciao world!
   Hello world!
   Bonjour world!

Installing Anaconda
-------------------

Pipeliner uses virtual environments managed by `Conda`, which is available through `Anaconda <https://www.continuum.io/downloads>`_. Download the distribution pre-packaged with Python 2.7. 

Make sure conda is installed and updated::

   conda --version
   conda update conda

.. tip:: If this is your first time working with `Conda`, you may need to edit your configuration paths to ensure `Anaconda` is invoked when calling ``conda``


Pre-Packaged Conda Environment
------------------------------

Yaml File
`````````

`Environment for Linux`::

   conda env create -f pipeliner/envs/linux_env.yml

`Environment for OS X`::

   conda env create -f pipeliner/envs/osx_env.yml


.. note:: Copies of pre-compiled binaries are hosted/maintained at https://anaconda.org/Pipeliner/repo

.. warning:: For those installing on the Shared Computing Cluster (SCC) at Boston University, instructions on how to setup a private conda environment can be `here <https://github.com/montilab/pipeliner/blob/master/scripts/mkenv.sh>`_.


Setting up Pipeliner
--------------------

With all prerequisites, one can quickly setup Pipeliner by cloning the repository, configuring local paths to toy datasets, activating the conda environment, and downloading the Nextflow executable::

   # Activate conda environment
   source activate pipeliner

   # Clone Pipeliner
   git clone https://github.com/montilab/pipeliner

   # Configure local paths to toy datasets
   python pipeliner/scripts/paths.py

   # Move to pipelines directory
   cd pipeliner/pipelines

   # Download nextflow executable
   curl -s https://get.nextflow.io | bash

   # Run RNA-seq pipeline with toy data
   ./nextflow rnaseq.nf -c rnaseq.config

The output should look like this::

   N E X T F L O W  ~  version 0.31.1
   Launching `rnaseq.nf` [nasty_pauling] - revision: cd3f572ab2
   [warm up] executor > local
   [31/1b2066] Submitted process > pre_fastqc (ggal_alpha)
   [23/de6d60] Submitted process > pre_fastqc (ggal_theta)
   [7c/28ee53] Submitted process > pre_fastqc (ggal_gamma)
   [97/9ad6c1] Submitted process > check_reads (ggal_alpha)
   [ab/c3eedf] Submitted process > check_reads (ggal_theta)
   [2d/050633] Submitted process > check_reads (ggal_gamma)
   [1d/f3af6d] Submitted process > pre_multiqc
   [32/b1db1d] Submitted process > hisat_indexing (genome_reference.fa)
   [3b/d93c6d] Submitted process > trim_galore (ggal_alpha)
   [9c/3fa50b] Submitted process > trim_galore (ggal_theta)
   [62/25fce0] Submitted process > trim_galore (ggal_gamma)
   [66/ccc9db] Submitted process > hisat_mapping (ggal_alpha)
   [28/69fff5] Submitted process > hisat_mapping (ggal_theta)
   [5c/5ed2b6] Submitted process > hisat_mapping (ggal_gamma)
   [b4/e559ab] Submitted process > gtftobed (genome_annotation.gtf)
   [bc/6f490c] Submitted process > rseqc (ggal_alpha)
   [71/80aa9e] Submitted process > rseqc (ggal_theta)
   [17/ca0d9f] Submitted process > rseqc (ggal_gamma)
   [d7/7d391b] Submitted process > counting (ggal_alpha)
   [df/936854] Submitted process > counting (ggal_theta)
   [11/143c2c] Submitted process > counting (ggal_gamma)
   [31/4c11f9] Submitted process > expression_matrix
   [1f/3af548] Submitted process > multiqc
   Success: Pipeline Completed!
