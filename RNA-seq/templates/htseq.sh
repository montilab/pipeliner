  samtools view ${bamfiles} | \\
  htseq-count \\
  --type '${params.htseq.type}' \\
  --mode '${params.htseq.mode}' \\
  --idattr '${params.htseq.idattr}' \\
  --order '${params.htseq.order}' \\
  - ${gtf} > '${sampleid}.counts.txt'