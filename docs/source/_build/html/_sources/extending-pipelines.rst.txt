.. _extending-pipelines-page:

Extending Pipelines
===================

.. toctree::
   :maxdepth: 3
   :caption: Contents:

General Workflow
----------------

The framework provides multiple resources for the user to extend and create sequencing pipelines. The first is toy datasets for all available pipelines including sequencing files, alignment files, genome reference and annotation files, as well as phenotypic data. Additionally, there are pre-defined scripts, processes, and configuration files that can be inherited and easily modified for various pipelines. Together, users can rapidly develop flexible and scalable pilelines. Lastly, there is a testing module enabling users to frequently test a series of different configurations with each change to the codebase.

Configuration Inheritance
-------------------------

An important property of configuration files is that they are inheritable. This allows developers to focus soley on the configuration components that are changing with each pipeline execution. Typically there are four components of a configuration file including the following.

Executor parameters::
    
    process {
        executor = "local"
    }

Input data file paths::

    indir  = "/Users/anthonyfederico/pipeliner/pipelines/toy_data/rna-seq"
    outdir = "/Users/anthonyfederico/pipeliner/pipelines/rna-seq-results"


Pipeline parameters::

    aligner      = "hisat"
    quantifier   = "htseq" 

Process-specific parameters::

    htseq.type   = "exon"
    htseq.mode   = "union"
    htseq.idattr = "gene_id"
    htseq.order  = "pos"

When developing, typically the only parameters that will be changing are pipeline parameters when testing the full scope of flexibility. Therefore, the development configuration file will look something like the following::

    // paired / hisat / featurecounts

    includeConfig "local.config"
    includeConfig "dataio.config"

    paired        = true
    aligner       = "hisat"
    quantifier    = "featurecounts"
    skip.counting = false
    skip.rseqc    = false
    skip.multiqc  = false
    skip.eset     = false

    includeConfig "parameters.config"

Template Process Injections
---------------------------

.. note:: Sometimes it's better to create a new template rather than heavily modify an existing one

Each pipeline is essentially a series of modules - connected through minimal Nextflow scripting - that execute pre-defined template processes. While templates are generally defined to be applicable to multiple pipelines and are parameterized in a configuration file, they have two additional components contributing to their flexibility.

The following is an example of a template process for the third-party software tool `featureCounts`:

.. code-block:: bash
   :linenos:

    featureCounts \\

    # Common flags directly defined by the user
    -T ${params.feature_counts.cpus} \\
    -t ${params.feature_counts.type} \\
    -g ${params.feature_counts.id} \\

    # Flags handled by the pipeline
    -a ${gtf} \\
    -o "counts.raw.txt" \\

    # Arguments indirectly defined by the user
    ${feature_counts_sargs} \\

    # Extra arguments
    ${params.feature_counts.xargs} \\

    # Input data
    ${bamfiles};

    # After injection 
    ${params.feature_counts.ainj}

**Lines 4-6**
  These are *common* keyword arguments that can be set to string/int/float types by the user and passed *directly* from the configuration file to the template. The *params* prefix in the variable means it is initialized in the configuration file.

**Lines 9-10**
  These are flags that are typically non-dynamic and handled interally by the pipeline.

**Line 13**
  These are *common* flags that must be *indirectly* defined by the user. For example, featurCounts requires a ``-p`` flag for paired reads. Because ``params.paired`` is a boolean, it makes more sense for the pipeline to create a string of supplemental arguments indirectly defined by the configuration file.

.. code-block:: bash

    feature_counts_sargs = ""
    if (params.paired) {
        feature_counts_sargs = feature_counts_sargs.concat("-p ")
    }

**Line 16**
  These are *uncommmon* keyword arguments or flags that can be pass *directly* from the configuration file to the template. Because some software tools can include hundreds of arguments, we explicitly state common arguments, but allow the user to additionally insert any unlimited number of additional arguments to maximize flexibility.

  For example, the user might want to perform a one-off test of the pipeline where they remove duplicate reads and only count fragments that have a length between 50-600 base pairs. These options can be injected into the template by simply defining ``params.feature_counts.xargs = "--ignoreDup -d 50 -D 600"`` in the configuration file.

**Line 19**
  These are required arguments such as input data handled interally by the pipeline.

**Line 22**
  These are code injections - typically one-liner cleanup commands - that can be injected after the main script of a template. For example, the output of featureCounts is a `genes x samples` matrix and the user may want to try sorting rows by gene names. Setting ``params.feature_counts.ainj`` to ``"sort -n -k1,1 counts.raw.txt > counts.raw.txt;"`` would accomplish such a task.


After parameterization, the final result would look something like this:

.. code-block:: bash
   :linenos:

    featureCounts -T 1 -t "exon" -g "gene_id" \
    -a "path/to/reference_annotation.gtf" \
    -o "counts.raw.txt" \
    -p --ignoreDup -d 50 -D 600 \
    s1.bam s2.bam s3 bam;
    sort -n -k1,1 counts.raw.txt > counts.raw.txt;

Testing Module
--------------
Each major change to a pipeline should be followed with a series of tests. Because pipelines are so flexible, it's infeasible to manually test even a limited set of typical configurations. To solve this problem we include an automated testing module.

Users can automatically test a series of configuration files by specifying a directory of user-defined tests::

    /pipeliner
     └── /tests
          └── /configs
               └── /rnaseq
                    ├── /t1.config
                    ├── /t2.config
                    └── /t3.config   
 

To run these series of tests, users can execute ``python pipeliner/test.py rnaseq`` which will search for the directory ``pipeliner/tests/configs/rnaseq`` and automatically pair and run each configuration file with a pipeline script named ``rnaseq.nf``.

.. note:: The directory name of tests must be the same as the pipeline script they are paired with
