hisat2 -x '${params.outdir}/hisat_files/hisat_index/part' -U ${reads[0]} \\
-S '${sampleid}.sam' --summary-file '${sampleid}.log' --new-summary
samtools view -S -b '${sampleid}.sam'
samtools sort '${sampleid}.sam' -o '${sampleid}.sorted.bam'