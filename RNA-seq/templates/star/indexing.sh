mkdir star_index
STAR --runMode genomeGenerate \\
--runThreadN ${params.star_indexing.cpus} \\
--sjdbGTFfile ${gtf} \\
--genomeDir star_index/ \\
--genomeFastaFiles ${fasta} \\
--sjdbOverhang ${params.star_indexing.sjdbOverhang}