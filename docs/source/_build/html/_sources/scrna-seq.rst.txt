.. _scrna-seq-page:

scRNA-seq
=========

.. toctree::
   :maxdepth: 3
   :caption: Contents:

Check Reads :code:`check_reads`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`)
:output: None
:script: Ensures correct format of sequencing read files

Genome Indexing :code:`hisat_indexing/star_indexing`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  Genome reference file (`.fa`) | Genome annotation file (`.gtf`)
:output: Directory containing indexed genome files
:script: Uses either `STAR` or `HISAT2` to build an indexed genome

Quality Check :code:`fastqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`)
:output: Report files (`.html`)
:script: Uses `FastQC` to check quality of read files

Whitelist :code:`whitelist`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`)
:output: A table of white listed barcodes (`.txt`)
:script: Uses `UMI-tools` to extract and identify true cell barcodes

Extract :code:`extract`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`) | A table of white listed barcodes (`.txt`)
:output: Extracted read files (`.fastq`)
:script: Uses `UMI-tools` to extract barcode from reads and append to read name

Read Mapping :code:`hisat_mapping/star_mapping`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`) | Genome annotation file (`.gtf`) | Directory containing indexed reference genome files
:output: A list of alignment files (`.bam`) | Log files (`.log`)
:script: Uses either `STAR` or `HISAT2` to align reads to a reference genome

Reformat Reference :code:`gtftobed`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: Genome annotation file (`.gtf`)
:output: Genome annotation file (`.bed`)
:script: Converts genome annotation file from GTF to BED format

Mapping Quality :code:`rseqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: A list of alignment files (`.bam`) 
:output: Report files (`.txt`)
:script: Uses `RSeQC` to check quality of alignment files

Quantification :code:`counting`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: A list of alignment files (`.bam`) | Genome annotation file (`.gtf`) 
:output: Read counts (`.txt`) | Log files (`.txt`)
:script: Uses `featureCounts` to quantify reads

Summary Report :code:`multiqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: Log files and summary reports from all processes
:output: A summary report (`.html`)
:script: Uses `MultiQC` to generate a summary report