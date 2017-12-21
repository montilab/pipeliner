mkdir index
STAR --runMode genomeGenerate \\
--runThreadN ${params.star_indexing.cpus} \\
--sjdbGTFfile ${gtf} \\
--genomeDir index/ \\
--genomeFastaFiles ${fasta} \\
--sjdbOverhang ${params.star_indexing.sjdb_overhang}