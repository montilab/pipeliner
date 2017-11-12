samtools index $bamfiles
bam_stat.py -i $bamfiles > ${sampleid}.bam_stats.txt
geneBody_coverage.py -r $bed -i $bamfiles -o ${sampleid} > ${sampleid}.read_coverage.txt
junction_annotation.py -i $bamfiles -o ${sampleid} -r $bed > ${sampleid}.read_junction.txt
read_distribution.py  -i $bamfiles -r $bed > ${sampleid}.read_distribution.txt
read_duplication.py -i $bamfiles -o ${sampleid}.read_duplication.txt