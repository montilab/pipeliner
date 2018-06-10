trim_galore \\
--gzip \\
--quality ${params.trim_galore.quality} \\
--fastqc \\
--adapter ${params.trim_galore.adapter1} \\
${reads[0]}