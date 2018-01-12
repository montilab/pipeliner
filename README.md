# Pipeliner   
<i>Flexible and robust framework for the specification of high-throughput sequencing data processing workflows</i>    

[![Built With](https://img.shields.io/badge/Built%20With-Nextflow-brightgreen.svg)](https://www.nextflow.io/)
![Python](https://img.shields.io/badge/Pipeline-Python%202.7-blue.svg)
![Python](https://img.shields.io/badge/Web%20App-Python%203.6-blue.svg)
![Compatibility](https://img.shields.io/badge/Compatibility-Linux%20%2F%20OSX-orange.svg)
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)
[![GitHub Issues](https://img.shields.io/github/issues/montilab/pipeliner.svg)](https://github.com/montilab/pipeliner/issues)

## Features   
* Modular directory structure: It is designed to generate automated result directory based on the names of the samples and tools used to process them   
* Platform independent: It is bundled with an anaconda repository which contains pre-compiled tools as well as pre-built environments that can use used directly.   
* Modular architecture: It allows the expert users to customize, modify processes, or add additional tools based on their needs.    
* Automated job parallelization, job recovery, and reproducibility.

## Quickstart
*For more information, please refer to the [full documentation](https://github.com/montilab/pipeliner/blob/master/docs/documentation.md)*

### Clone Repository
```bash
$ git clone https://github.com/montilab/pipeliner
```

### Create Conda Environment and Install Dependencies
```bash
$ conda env create -f pipeliner/RNA-seq/envs/linux_env.yml
```

```bash
$ conda env create -f pipeliner/RNA-seq/envs/osx_env.yml
```

### Activate Conda Environment
```bash
$ source activate pipeliner
```

### Update Local Paths
```bash
$ python pipeliner/scripts/paths.py
```

### Download Nextflow Executable
```
$ cd pipeliner/RNA-seq
$ curl -s https://get.nextflow.io | bash
```

### Locally Run Example Data
```bash
$ ./nextflow main.nf -c nextflow.config
```

### Expected Output
```text
Launching `main.nf` [distraught_hugle] - revision: 0e9a7a8940
 P I P E L I N E R  ~  v2.4
====================================
Reads          : /Users/anthonyfederico/pipeliner/RNA-seq/ggal_data/ggal_reads.csv
Reference      : /Users/anthonyfederico/pipeliner/RNA-seq/ggal_data/genome_reference.fa
Annotation     : /Users/anthonyfederico/pipeliner/RNA-seq/ggal_data/genome_annotation.gtf
Input Dir      : /Users/anthonyfederico/pipeliner/RNA-seq/ggal_data
Output Dir     : /Users/anthonyfederico/pipeliner/RNA-seq/ggal_results
====================================
Read Type      : paired-end
Aligner        : hisat
Quantifier     : htseq
Save Reference : true
Save Temporary : true
====================================
Current user  : anthonyfederico
Current home  : /Users/anthonyfederico
Current path  : /Users/anthonyfederico/pipeliner/RNA-seq
====================================
[warm up] executor > local
[32/b1db1d] Submitted process > hisat_indexing (genome_reference.fa)
[3b/d93c6d] Submitted process > trim_galore (ggal_alpha)
[9c/3fa50b] Submitted process > trim_galore (ggal_theta)
[62/25fce0] Submitted process > trim_galore (ggal_gamma)
[96/2ffc07] Submitted process > fastqc (ggal_alpha)
[31/89cddb] Submitted process > fastqc (ggal_theta)
[57/b2fdf0] Submitted process > fastqc (ggal_gamma)
[66/ccc9db] Submitted process > hisat_mapping (ggal_alpha)
[28/69fff5] Submitted process > hisat_mapping (ggal_theta)
[5c/5ed2b6] Submitted process > hisat_mapping (ggal_gamma)
[bc/6f490c] Submitted process > rseqc (ggal_alpha)
[71/80aa9e] Submitted process > rseqc (ggal_theta)
[17/ca0d9f] Submitted process > rseqc (ggal_gamma)
[d7/7d391b] Submitted process > htseq (ggal_alpha)
[df/936854] Submitted process > htseq (ggal_theta)
[11/143c2c] Submitted process > htseq (ggal_gamma)
[1f/3af548] Submitted process > multiqc
Success: Pipeline Completed!
```

*Please refer to [reports](https://github.com/montilab/pipeliner/blob/master/docs/reports.md) for examples.*