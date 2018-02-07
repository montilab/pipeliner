stringtie ${bamfiles} \\
-C '${sampleid}.covrefs.gtf' \\
-o '${sampleid}.transcripts.gtf' \\
-A '${sampleid}.counts.txt' \\
-v -G ${gtf}