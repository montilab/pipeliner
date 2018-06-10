mkdir index;
hisat2-build -p ${params.hisat_indexing.cpus} \\
             -f $fasta \\
             index/part