trim_galore \\
${params.trim_galore.xargs} \\
--gzip \\
--quality ${params.trim_galore.quality} \\
--fastqc \\
--adapter ${params.trim_galore.adapter1} --adapter2 ${params.trim_galore.adapter2} \\
--paired ${reads[0]} ${reads[1]}