bowtie2 -x bowtie_index -1 ${reads[0]} -2 ${reads[1]} -S output.sam
samtools view -Sb output.sam > output.bam
samtools sort output.bam -o output.sorted.bam