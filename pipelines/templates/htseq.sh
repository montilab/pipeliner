samtools view ${bamfiles} | \\
htseq-count \\
${params.htseq.xargs} \\
--type '${params.htseq.type}' \\
--mode '${params.htseq.mode}' \\
--idattr '${params.htseq.idattr}' \\
--order '${params.htseq.order}' \\
- ${gtf} > '${sampleid}.counts.txt';
${params.htseq.ainj}