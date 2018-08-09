.. _dege-page:

DGE
===

.. toctree::
   :maxdepth: 3
   :caption: Contents:

Quantification :code:`counting`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: A list of alignment files (`.bam`) | Genome annotation file (`.gtf`) 
:output: Read counts (`.txt`) | Log files (`.txt`)
:script: Uses `featureCounts` to quantify reads

Expression Matrix :code:`expression_matrix`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: A list of count files (`.txt`)
:output: An expression matrix (`.txt`)
:script: Reformats a list of count files into a genes x samples matrix

Sample Renaming :code:`rename_samples`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: An expression matrix (`.txt`)
:output: An expression matrix (`.txt`)
:script: Renames samples in expression matrix based on a user-supplied table

Summary Report :code:`multiqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: Log files and summary reports from all processes
:output: A summary report (`.html`)
:script: Uses `MultiQC` to generate a summary report