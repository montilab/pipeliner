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
  Markdown (MD) and Restructured Text (reST) documentaion files associated with Pipeliner and existing pipelines

**envs**
  Files (YAML) and scripts requires to reproduce Conda environments

**scripts**
  Various helper scripts for framework setup and maintainence

**tests**
  Testing module (Python) for multi-pipeline automatic test execution and reporting

**pipelines/configs**
  Base config files inherited by pipeline configurations

**pipelines/scripts**
  Various helper scripts for pipeline processes

**pipelines/templates**
  Template processes inherited by pipeline workflows

**pipelines/toy_data**
  Small datasets for rapid development and testing of pipelines

**pipelines/rnaseq.nf**
  Nextflow script for the RNA-seq pipeline

**pipelines/rnaseq.config**
  Configuration file for the RNA-seq pipeline
