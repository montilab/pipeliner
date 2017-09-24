wget ftp://ftp.ensembl.org/pub/release-79/gtf/homo_sapiens/Homo_sapiens.GRCh38.79.gtf.gz
gunzip Homo_sapiens.GRCh38.79.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-79/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://www.encodeproject.org/files/ENCFF001RFH/@@download/ENCFF001RFH.fastq.gz -O ENCFF001RFH.fastq.gz 
gunzip ENCFF001RFH.fastq.gz 

wget https://www.encodeproject.org/files/ENCFF001RFG/@@download/ENCFF001RFG.fastq.gz -O ENCFF001RFG.fastq.gz
gunzip ENCFF001RFG.fastq.gz
