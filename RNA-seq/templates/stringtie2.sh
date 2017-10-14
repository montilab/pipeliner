stringtie ${bamfiles} \\
-o ${bamfiles}_transcripts.gtf \\
-v \\
-G ${mergedgtf} \\
-A ${bamfiles}.gene_abund.txt \\
-C ${bamfiles}.cov_refs.gtf \\
-e \\
-b ${bamfiles}_ballgown