samtools index $bamfiles
bam_stat.py -i $bamfiles > ${sampleid}.bam_stats
geneBody_coverage.py -r $bed -i $bamfiles -o ${sampleid}
junction_annotation.py -i $bamfiles -o ${sampleid} -r $bed