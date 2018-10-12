trim_galore \\
${params.trim_galore.xargs} \\
--gzip \\
--quality ${params.trim_galore.quality} \\
--fastqc \\
${reads[0]}