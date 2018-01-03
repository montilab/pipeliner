samtools index $bamfiles
bam_stat.py -i $bamfiles > ${sampleid}.bam_stats