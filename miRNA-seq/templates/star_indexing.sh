mkdir star_index
STAR --runMode genomeGenerate \\
--runThreadN ${task.cpus} \\
--sjdbGTFfile ${gtf} \\
--genomeDir star_index/ \\
--genomeFastaFiles ${fasta} \\
--sjdbOverhang ${params.star_indexing.sjdbOverhang}