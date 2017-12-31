featurecounts \\
-p \\
-t '${params.feature_counts.type}' \\
-g '${params.feature_counts.id}' \\
-a ${gtf} -o '${sampleid}.counts.txt' ${bamfiles}