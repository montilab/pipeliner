samtools index $bamfiles
bam_stat.py -i $bamfiles > ${sampleid}.bam_stats.summary.txt
geneBody_coverage.py -r $bed -i $bamfiles -o ${sampleid} > ${sampleid}.read_coverage.summary.txt
junction_annotation.py -i $bamfiles -o ${sampleid} -r $bed > ${sampleid}.read_junction.summary.txt
read_distribution.py  -i $bamfiles -r $bed > ${sampleid}.read_distribution.summary.txt
read_duplication.py -i $bamfiles -o ${sampleid}.read_duplication.summary.txt