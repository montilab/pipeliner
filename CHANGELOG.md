### Pipeliner (RNA-Seq)
### *Change Log*

### Version 2.0

- add introspection to parameter log
	- timestamp
	- work directory
	- command
- create a more informative console message on startup
- update nextflow.config with standard and process-specific paramaters
- update nextflow.config with comments
- separated out part of rseqc module temporarily
- introduced bool `from_bam` to start directly from bamfiles
- aggrcounts.py -> aggregator.py

### Version 1.1

- introduced bowtie2 as an alignment option for simple organisms
- changed variable `star_index` to `index` to make the boolean generic to any aligner
- separated out kallisto module temporarily
- created `afterProcess{}` for optional temp file removal
- changed deprecated `.first` and `.into` channel operators
- updated documentation for working with the shared computing cluster 
- added full tutorial for simple RNA-seq pipeline with toy dataset
- updated documentation for setting up nexflow
- separate conda environments .yml created for linux and mac
- switched from miniconda to anaconda
- separated custom python files into a scripts folder
- templates are now dynamically generated using key word arguements set in settings.py
- moved all module calls to nextflow script templates
- web app runs Python 3.6 while any environment can be used for nextflow
- flask sockets enabled for two way communication between web console and nextflow
- web app dynamically looks up existing files and directories
- flask api now wraps around nextflow
- reduced duplicate variables for aligning index
- created prototype of web application made with flask framework and semantic styling
- `zcat` changed to `gunzip -c` for linux/mac support 
- simplified configuration file structure
- changed flowing data variable names to be more orderly
- removed *most* hard coded file paths
	- infile/outfile remain hardcoded
	- fastq reads in .csv file remain hardcoded

### Version 1.0