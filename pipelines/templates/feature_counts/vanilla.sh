featureCounts \\
$feature_counts_sargs \\
${params.feature_counts.xargs} \\
-T ${params.feature_counts.cpus} \\
-t ${params.feature_counts.type} \\
-g ${params.feature_counts.id} \\
-a ${gtf} \\
-o 'counts.raw.txt' \\
${bamfiles};
${params.feature_counts.ainj}