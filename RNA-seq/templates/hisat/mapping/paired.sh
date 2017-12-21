hisat2 -p ${params.hisat_mapping.cpus} --summary-file '${sampleid}.hisat.log' --new-summary -x '${params.outdir}/hisat/index/part' -1 ${reads[0]} -2 ${reads[1]} -S '${sampleid}.sam'
samtools view -S -b '${sampleid}.sam'
samtools sort '${sampleid}.sam' -o '${sampleid}.sorted.bam'