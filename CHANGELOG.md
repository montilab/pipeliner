### Pipeliner (RNA-Seq)
### *Change Log*

### Version 1.1

- separated out kallisto temporarily
- created `afterProcess{}` for optional temp file removal
- changed deprecated `.first` and `.into` channel operators
- updated documentation for working with the shared computing cluster 
- full tutorial for simple RNA-seq pipeline with toy dataset
- updated documentation for setting up nexflow
- separate conda environments .yml created for linux and mac
- switch from miniconda to anaconda
- separated custom python files into a scripts folder
- templates are dynamically generated using key word arguements set in settings.py
- moved all module calls to nextflow script templates
- webapp runs Python 3.6 while any environment can be used for nextflow
- flask sockets enabled for two way communication between web console and nextflow
- dynamic file lookups within web app
- flask api wraps around nextflow
- reduced duplicate variables for aligning index
- prototype of web application made with flask and semantic
- `zcat` changed to `gunzip -c` for linux/mac support 
- simplified configuration file structure
- changed flowing data variable names to be more orderly
- removed *most* hard coded file paths
	- infile/outfile remain hardcoded
	- fastq reads in .csv file remain hardcoded

### Version 1.0