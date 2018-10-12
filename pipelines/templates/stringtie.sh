stringtie ${bamfiles} \\
${params.stringtie.xargs} \\
-C '${sampleid}.covrefs.gtf' \\
-o '${sampleid}.transcripts.gtf' \\
-A '${sampleid}.counts.txt' \\
-v -G ${gtf};
${params.stringtie.ainj}