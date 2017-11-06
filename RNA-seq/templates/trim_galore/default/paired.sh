trim_galore \\
--gzip \\
--quality ${params.trim_galore.quality} \\
--paired ${reads[0]} ${reads[1]}