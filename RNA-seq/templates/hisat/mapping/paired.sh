hisat2 -x '${params.outdir}/hisat_files/hisat_index/part' -1 ${reads[0]} -2 ${reads[1]} \\
-S '${sampleid}.sam' --summary-file '${sampleid}.log' --new-summary
samtools view -S -b '${sampleid}.sam'
samtools sort '${sampleid}.sam' -o '${sampleid}.sorted.bam'