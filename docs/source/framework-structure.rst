.. _framework-structure-page:

Framework Structure
===================

.. toctree::
   :maxdepth: 3
   :caption: Contents:

Pipeline is a framework with various moving parts to support the development of multiple sequencing pipelines. The following is a simplified example of its directory structure::

    /pipeliner
     ├── /docs
     ├── /envs
     ├── /scripts
     ├── /tests
     └── /pipelines
          ├── /configs
          ├── /scripts
          ├── /templates
          ├── /toy_data
          ├── /rnaseq.nf
          └── /rnaseq.config     

**docs**
  Markdown and Restructured Text documentaion files associated with Pipeliner and existing pipelines

**envs**
  Yaml files and scripts required to reproduce Conda environments

**scripts**
  Various helper scripts for framework setup and maintenance

**tests**
  Python testing module for multi-pipeline automatic test execution and reporting

**pipelines/configs**
  Base config files inherited by pipeline configurations

**pipelines/scripts**
  Various helper scripts for pipeline processes

**pipelines/templates**
  Template processes inherited by pipeline workflows

**pipelines/toy_data**
  Small datasets for rapid development and testing of pipelines. These datasets are modifications from original `RNA-seq <https://github.com/nextflow-io/rnatoy/tree/master/data/ggal>`_ and `scRNA-seq <http://cf.10xgenomics.com/samples/cell-exp/1.3.0/hgmm_100/hgmm_100_fastqs.tar>`_ datasets.

**pipelines/rnaseq.nf**
  Nextflow script for the RNA-seq pipeline

**pipelines/rnaseq.config**
  Configuration file for the RNA-seq pipeline
