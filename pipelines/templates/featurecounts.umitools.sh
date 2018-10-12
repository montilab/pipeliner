featureCounts \\
${params.feature_counts.xargs} \\
-a ${gtf} \\
-o '${sampleid}.counts.txt' \\
-R BAM ${bamfiles} \\
-T 4;  
samtools sort '${sampleid}.bam.featureCounts.bam' -o assigned_sorted.bam;
samtools index assigned_sorted.bam;
umi_tools count --per-gene --gene-tag=XT --per-cell -I assigned_sorted.bam -S '${sampleid}.counts.tsv.gz';
${params.feature_counts.ainj}