umi_tools extract --bc-pattern=${params.extract.bc_pattern} \\
--stdin ${reads[0]} \\
--stdout ${sampleid}_R1_extracted.fastq.gz \\
--read2-in ${reads[1]} \\
--read2-out=${sampleid}_R2_extracted.fastq.gz \\
--filter-cell-barcode \\
--whitelist='${params.outdir}/${sampleid}/whitelisted/whitelist.txt';