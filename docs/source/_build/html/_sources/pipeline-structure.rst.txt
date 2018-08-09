.. _pipeline-structure-page:

Pipeline Structure
==================

.. toctree::
   :maxdepth: 3
   :caption: Contents:

Input
-----
The file paths for all data fed to a pipeline are specified in the configuration file. To ease the development process, Pipeline provides toy datasets for each of the pipelines. This example will cover the RNA-seq pipeline.


.. note:: Data for this pipeline is located in **pipelines/toy_data/rna-seq**

Users must provide the following files:
    - Sequencing files or alignment files
    - Comma-delimited file containing file paths to reads/bams
    - Genome reference file
    - Genome annotation file

Configuration File
------------------

The configuration file is where all file paths are specified and pipeline processes are paramaterized. The configuration can be broken into three sections, including file paths, executor and compute resources, and pipeline options and parameters.

File Paths
``````````
The configuration file specifies where to find all of the input data. Additionally, it provide a path to an output directory where the pipeline will output results. The following is a typical example for the RNA-seq configuration file::

    indir  = "/Users/anthonyfederico/pipeliner/pipelines/toy_data/rna-seq"
    outdir = "/Users/anthonyfederico/pipeliner/pipelines/rna-seq-results"
    fasta  = "${params.indir}/genome_reference.fa"
    gtf    = "${params.indir}/genome_annotation.gtf"
    reads  = "${params.indir}/ggal_reads.csv"

Executor and Compute Resources
``````````````````````````````

An abstraction layer between Nextflow and Pipeliner logic enables platform independence and seamless compatibility with high performance computing executors. This allows users to execute pipelines on their local machine or through a computing cluster by simply specifying in the configuration file.

Pipeliner provides two base configuration files that can be inherited depending if a pipeline is being executing using local resources or a Sun Grid Engine (SGE) queuing system. 

If the latter is chosen, pipeline processes will be automatically parallelized. Additionally, each individual process can be allocated specific computing resource instructions when nodes are requested.

*Local config example*::

    process {
      executor = 'local'
    }

*Cloud computing config example*::

    process {
      executor = 'sge'
      scratch = true

      $trim_galore.clusterOptions            = "-P montilab -l h_rt=24:00:00 -pe omp 8"
      $star_mapping.clusterOptions           = "-P montilab -l h_rt=24:00:00 -l mem_total=94G -pe omp 16"
      $counting.clusterOptions               = "-P montilab -l h_rt=24:00:00 -pe omp 8"
      $expression_matrix.clusterOptions      = "-P montilab -l h_rt=24:00:00 -pe omp 8"
      $multiqc.clusterOptions                = "-P montilab -l h_rt=24:00:00 -pe omp 8"
    }

Pipeline Options and Parameters
```````````````````````````````

The rest of the configuration file is dedicated to the different pipeline options and process parameters that can be specified. Some import examples include the following::

    # General pipeline parameters
    aligner      = "hisat"
    quantifier   = "htseq" 

    # Process-specific parameters
    htseq.type   = "exon"
    htseq.mode   = "union"
    htseq.idattr = "gene_id"
    htseq.order  = "pos"    

Pipeline Script
---------------

Template Processes
------------------

Pipelines written in Nextflow consist of a series of processes. Processes specify data I/O and typically wrap around third-party software tools to process this data. Processes are connected through channels – asynchronous FIFO queues – which manage the flow of data throughout the pipeline.

Processes have the following basic structure::
    
    process <name> {

        input:
        <process inputs>

        output:
        <process outputs>

        script:
        <user script to be executed>
    }


Often, the script portion of the processes are reused by various sequencing pipelines. To help standardize pipeline development and ensure good practices are propogated to all pipelines, template processes are defined and inherited by pipeline processes.

.. note:: Templates are located in **pipelines/templates**

For example, these two processes execute the same code::

    # Without inheritance
    process htseq {

        input:
        <process inputs>

        output:
        <process outputs>

        script:
        '''
        samtools view ${bamfiles} | htseq-count - ${gtf} \\
        --type   ${params.htseq.type}   \\
        --mode   ${params.htseq.mode}   \\
        --idattr ${params.htseq.idattr} \\
        --order  ${params.htseq.order}  \\
        > counts.txt
        '''
    }

    # With inheritance
    process htseq {

        input:
        <process inputs>

        output:
        <process outputs>

        script:
        template 'htseq.sh'
    }

Output
------

The RNA-seq pipeline output has the following basic structure::

    /pipeliner/RNA-seq
    └── /results
        │
        ├── /sample_1
        │   ├── /trimgalore      | Trimmed Reads (.fq.gz) for sample_1
        │   ├── /fastqc 
        │   ├── /rseqc          
        │   └── /htseq
        │
        ├── /alignments          | Where (.bams) are saved
        ├── /aligner             
        │   └── /index           | Index created and used during mapping
        │
        ├── /expression_matrix   | Aggregated count matrix
        ├── /expression_set      | An expression set (.rds) object
        ├── /reports             | Aggregated report across all samples pre/post pipeliner
        └── /logs                | Process-related logs


Each sample will have its own directory with sample-specific data and results for each process. Additionally, sequencing alignment files and the indexed reference genome will be saved for future use if specified. Summary reports pre/post-workflow can be found inside the reports directory.
