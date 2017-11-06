mkdir hisat_index
hisat2-build -p ${params.hisat_indexing.cpus} -f $fasta hisat_index/part