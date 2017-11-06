# Pre-requisites
Pipeliner requires Java, Nextflow, and Anaconda for implementation. All other tools for implementation are wrapped in the conda environment described below. 

#### Test Nextflow

Make sure you have Java 7/8 installed and then install Nextflow to any working directory. Test the Nextflow executable before continuuing. This is just to make sure your environment is compatible with a nextflow executable. You will download another one later in the tutorial.
```bash
java -version
cd path/to/wd
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

#### Download Conda

*Local Machine*
Conda is available through [Anaconda](https://www.continuum.io/downloads). Download the distribution pre-packaged with Python 2.7. If this is your first time working with conda, you may need to edit your configuration paths to ensure anaconda is invoked when calling `conda`.

*Shared Computing Cluster*
Enable conda by loading a pre-installed version of Anaconda with `module load anaconda2`. This will load the latest pre-installed version pre-packaged with Python 2.7.

#### Create Conda Environment

```bash
conda env create -f envs/linux_env.yml 
source activate pipeliner

or

conda env create -f envs/osx_env.yml 
source activate pipeliner
```

# Setting Up Pipeliner

#### Clone Repository and Download Nextflow Executable

```bash
git clone https://github.com/montilab/pipeliner
cd pipeliner/RNA-seq
curl -s https://get.nextflow.io | bash
```
#### Pipeliner Structure
```text
/pipeliner/RNA-seq
├── /data
│   │
│   ├── /alignments
│   │   └── sample_1.bam       [^]
│   ├── /reads
│   │   ├── sample_1_R1.fq.gz  [^]
│   │   └── sample_1_R2.fq.gz  [^]
│   │
│   ├── genome_annotation.gtf  [^]
│   ├── genome_reference.fa    [^]
│   ├── genome_refseq.bed      [^]
│   ├── reads.csv              [*]
│   └── alignments.csv         [*]
│
├── nextflow
├── main.nf                    [*]
├── local.config               [*]
├── cluster.config             [*]
│
├── /templates
└── /scripts
```
[^] - Files that must be **uploaded** by the user  
[*] - Files that must be **modified** by the user

#### Files to Upload
```text
/data/reads
/data/alignments
/data/genome_annotation.gtf
/data/genome_reference.fa
/data/genome_refseq.fa 
```

#### Files to Modify
```text
/data/reads.csv
/data/alignments.csv
main.nf
local.config
cluster.config
```

#### Other Files and Folders
```text
nextflow
/templates
/scripts
```

# Running Pipeliner

#### Run
```bash
./nextflow main.nf -c local.config
```

#### Expected Output
```text
Launching `main.nf` [distraught_hugle] - revision: 0e9a7a8940
 P I P E L I N E R  ~  v2.3
====================================
Reads          : /Users/anthonyfederico/Village/pipeliner/RNA-seq/ggal_data/ggal_reads.csv
FASTA          : /Users/anthonyfederico/Village/pipeliner/RNA-seq/ggal_data/genome_reference.fa
Annotation     : /Users/anthonyfederico/Village/pipeliner/RNA-seq/ggal_data/genome_annotation.gtf
Input Dir      : /Users/anthonyfederico/Village/pipeliner/RNA-seq/ggal_data
Output Dir     : /Users/anthonyfederico/Village/pipeliner/RNA-seq/ggal_results
====================================
Read Type      : paired-end
Aligner        : star
Quantifier     : stringtie
Save Reference : true
Save Temporary : true
====================================
Current user  : anthonyfederico
Current home  : /Users/anthonyfederico
Current path  : /Users/anthonyfederico/Village/pipeliner/RNA-seq
====================================
[warm up] executor > local
[2d/f9d761] Submitted process > star_indexing (genome_reference.fa)
[b5/850c51] Submitted process > trim_galore (ggal_alpha)
[55/013bb2] Submitted process > fastqc (ggal_alpha)
[78/cad1fe] Submitted process > star_mapping (ggal_alpha)
[ed/72ccfa] Submitted process > stringetie (ggal_alpha)
[d5/e8494d] Submitted process > rseqc (ggal_alpha)
[f0/8653fb] Submitted process > stringtie_merge
[b0/b688b5] Submitted process > stringtie_eb (ggal_alpha)
[52/dbaf09] Submitted process > aggregate_counts
[55/9e3e61] Submitted process > multiqc
Success: Pipeline Completed!
```

#### Resume Feature
If there is an error with your workflow, you can fix it and return where you left off with
```
./nextflow main.nf -resume
```
