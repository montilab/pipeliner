mkdir star_index
STAR --runMode genomeGenerate --runThreadN ${task.cpus} --sjdbGTFfile $gtf --genomeDir star_index/ --genomeFastaFiles $fasta --sjdbOverhang 149