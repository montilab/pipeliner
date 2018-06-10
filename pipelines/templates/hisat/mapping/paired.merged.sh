hisat2 -p ${params.hisat_mapping.cpus} \\
     --summary-file '${sampleid}.log' \\
     --new-summary \\
     -x $index_base \\
     -U ${reads[1]} \\
     -S '${sampleid}.sam';
samtools sort '${sampleid}.sam' -o '${sampleid}.bam';