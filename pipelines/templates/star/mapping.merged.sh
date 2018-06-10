STAR --genomeDir $index \\
--sjdbGTFfile $gtf \\
--readFilesIn ${reads[1]} \\
--runThreadN ${params.star_mapping.cpus} \\
--readFilesCommand gunzip -c \\
--outFilterType BySJout \\
--outSAMtype BAM SortedByCoordinate \\
--outWigType bedGraph \\
--twopassMode ${params.star_mapping.twopassMode} \\
--outFilterMultimapNmax ${params.star_mapping.outfilter_multimap_nmax} \\
--outFilterMismatchNmax ${params.star_mapping.outfilter_mismatch_nmax} \\
--outFilterMismatchNoverLmax ${params.star_mapping.outfilter_mismatch_relmax} \\
--alignIntronMin ${params.star_mapping.align_intron_min} \\
--alignIntronMax ${params.star_mapping.align_intron_max} \\
--alignMatesGapMax ${params.star_mapping.align_mates_gapmax} \\
--alignSJoverhangMin ${params.star_mapping.align_sjoverhang_min} \\
--outFileNamePrefix "$sampleid.";
mv '${sampleid}.Aligned.sortedByCoord.out.bam' '${sampleid}.bam';