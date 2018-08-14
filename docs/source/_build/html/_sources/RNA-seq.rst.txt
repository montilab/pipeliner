.. _rna-seq-page:

RNA-seq
=======

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

Pre-Quality Check :code:`pre_fastqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`)
:output: Report files (`.html`)
:script: Uses `FastQC` to check quality of read files

Pre-MultiQC :code:`pre_multiqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  Log files (`.log`)
:output: Summary report file (`.html`)
:script: Uses `MultiQC` to generate a summary report

Read Trimming :code:`trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input:  List of read files (`.fastq`)
:output: Trimmed read files (`.fastq`) | Report files (`.html`)
:script: Trims low quality reads with `TrimGalore` and checks quality with `FastQC`

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
:script: Uses either `StringTie`, `HTSeQ`, or `featureCounts` to quantify reads

Expression Matrix :code:`expression_matrix`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: A list of count files (`.txt`)
:output: An expression matrix (`.txt`)
:script: Reformats a list of count files into a `genes x samples` matrix

Expression Features :code:`expression_features`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: Genome annotation file (`.gtf`) | An expression matrix (`.txt`)
:output: Gene feature data (`.txt`)
:script: Parses the genome annotation file for gene feature data

Expression Set :code:`expression_set`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: An expression matrix (`.txt`) | Gene feature data (`.txt`) | Sample phenotypic data (`.txt`)
:output: An expression set object (`.rds`)
:script: Creates an expression set object with eData, fData, and pData attributes

Summary Report :code:`multiqc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:input: Log files and summary reports from all processes
:output: A summary report (`.html`)
:script: Uses `MultiQC` to generate a summary report