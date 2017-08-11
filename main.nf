#!/usr/bin/env nextflow

// Pipeline version
version = 1.1

// Configurable variables
params.indir  = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data'
params.outdir = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_results'

params.fasta = "${params.indir}/genome_reference.fa"
params.gtf   = "${params.indir}/genome_annotation.gff"
params.reads = "${params.indir}/ggal_alpha.csv"
params.bed   = "${params.indir}/ggal_refseq.bed"

params.aligner = 'star'
params.star_index = false
params.paired = true
params.save_reference = true


// ---------------------------------------------//
//                VALIDATE FILES                //
// ---------------------------------------------//

// Read and map reads with samples using csv file
if (params.paired) {
  Channel
    .fromPath(params.reads)
    .splitCsv(header: true)
    .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
    .ifEmpty {error "File ${params.reads} not parsed properly"}
    .into {reads_fastqc; reads_trimgalore}
} else {
  Channel
    .fromPath(params.reads)
    .splitCsv(header: true)
    .map {row -> [row.Sample_Name, [row.Read]]}
    .ifEmpty {error "File ${params.reads} not parsed properly"}
    .into {reads_fastqc; reads_trimgalore}
}

// Validate inputs
if (params.star_index && params.aligner == 'star') {
  star_index = Channel
    .fromPath(params.star_index)
    .ifEmpty {exit 1, "STAR index not found: ${params.star_index}"}
    .toList()
}
if (params.gtf) {
  Channel
    .fromPath(params.gtf)
    .ifEmpty {exit 1, "GTF annotation file not found: ${params.gtf}"}
    .toList()
    .into {gtf_indexing; gtf_mapping; gtf_stringtie;}
}
if (params.bed) {
  bed = Channel
    .fromPath(params.bed)
    .ifEmpty {exit 1, "BED annotation file not found: ${params.bed}"}
    .toList()
}



log.info " P I P E L I N E R  ~  v${version}"
log.info "===================================="
log.info "Reads         : ${params.reads}"
log.info "Genome        : ${params.genome}"
log.info "FASTA         : ${params.fasta}"
log.info "Annotation    : ${params.gtf}"
log.info "Input Dir     : ${params.indir}"
log.info "Output Dir    : ${params.outdir}"
log.info "Aligner       : ${params.aligner}"
if (params.star_index) {
  log.info "Index         : ${params.star_index}"
}
log.info "Current user  : $USER"
log.info "Current home  : $HOME"
log.info "Current path  : $PWD"
log.info "===================================="

log_params(params, "nextflow_paramters.txt")

def log_params(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"
  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

// ---------------------------------------------//
//                BEGIN PIPELINE                //
// ---------------------------------------------//

// FASTQC 
process fastqc {
  tag "$sampleid"
  cache 'deep'
  publishDir "${params.outdir}/$sampleid/fastqc", mode: 'copy'

  input:
  set sampleid, reads from reads_fastqc

  output:
  set sampleid, '*_fastqc.{zip,html}' into fastqc_results

  script:
  if (params.paired) {
    """
    fastqc -o . -q ${reads[0]} ${reads[1]}
    """
  } else {
    """
    fastqc -o . -q ${reads[0]}
    """
  }
}

// Delete reads_fastqc

// TRIM GALORE
process trim_galore {
  tag "$sampleid"
  cache 'deep'
  publishDir "${params.outdir}/$sampleid/trimgalore", mode: 'copy'

  input:
  set sampleid, reads from reads_trimgalore

  output:
  set sampleid, '*fq.gz' into trimmed_reads
  set sampleid, '*trimming_report.txt' into trimgalore_results

  script:
  if(params.paired)
  {
    """
    trim_galore --paired --gzip ${reads[0]} ${reads[1]}
    """
  } else {
    """
    trim_galore --gzip ${reads[0]}
    """
  }
}

// Delete reads_trimgalore

// STAR - BUILD INDEX
if (params.aligner == 'star' && !params.star_index && file(params.fasta)) {
  process star_indexing {

    publishDir path: "${params.outdir}/star_files", saveAs: { params.save_reference ? it : null }, mode: 'copy'
    tag "$fasta"
    cache 'deep'

    input:
    file fasta from file(params.fasta)
    file gtf from gtf_indexing

    output:
    file "star_index" into star_index

    script:
    """
    mkdir star_index
    STAR \\
    --runMode genomeGenerate \\
    --runThreadN ${task.cpus} \\
    --sjdbGTFfile $gtf \\
    --sjdbOverhang 149 \\
    --genomeDir star_index/ \\
    --genomeFastaFiles $fasta
    """
  }
}

// STAR - MAP READS
if (params.aligner == 'star') {
  process star_mapping {
    tag "$sampleid"
    cache 'deep'
    publishDir "${params.outdir}/$sampleid/star", mode: 'copy'

    input:
    file index from star_index.first()
    file gtf from gtf_mapping.first()
    set sampleid, file (reads:'*') from trimmed_reads

    output:
    set sampleid, file("*.bam") into bam_files, bam_counts, bam_stringtie1, bam_stringtie2, bam_rseqc
    set sampleid, '*.out' into alignment_logs
    set sampleid, '*.tab' into alignment_tab

    script:
    """
    STAR \\
    --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --readFilesIn $reads  \\
    --runThreadN ${task.cpus} \\
    --twopassMode Basic \\
    --outWigType bedGraph \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand gunzip -c \\
    --outFileNamePrefix \'${sampleid}.'
    """
  }
}

// RSEQC
def num_bams
bam_counts.count().subscribe{num_bams=it}

process rseqc {
  tag "$sampleid"
  publishDir "${params.outdir}/$sampleid/rseqc", mode: 'copy'

  input:
  set sampleid, file(bamfiles) from bam_rseqc
  file bed from file(params.bed)

  output:
  file ("*.{txt,pdf,r,xls}") into rseqc_results
  file ('*.bam_stats') into bam_stats_results
  file ('*geneBodyCoverage.*') into gene_coverage_results
  stdout into gene_coverage_log
  file ('*junction.*') into junction_annotation_results
  //file('junc_annot.junction.xls') into junction
  //file('gene_coverage.geneBodyCoverage.txt') into coverage
    
  script:
  """
  samtools index $bamfiles
  bam_stat.py -i $bamfiles > ${sampleid}.bam_stats
  geneBody_coverage.py -r $bed -i $bamfiles -o ${sampleid}
  junction_annotation.py -i $bamfiles -o ${sampleid} -r $bed
  """
}

// STRINGTIE
if(params.aligner == "star" | params.aligner == "bowtie") {
  process stringtie1 {
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/stringtie1", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_stringtie1
    file gtf from file(params.gtf)

    output:
    file '*_transcripts.gtf' into gtf_list

    script:
    """
    stringtie $bamfiles \\
    -o ${bamfiles}_transcripts.gtf \\
    -v \\
    -G $gtf \\
    """
  }
  process gtf_mering {
    publishDir "${params.outdir}/stringtiemerge", mode: 'copy'

    input:
    val gtfs from gtf_list.toList()
    file gtf from file(params.gtf)

    output:
    file 'merged.gtf' into merged_gtf

    script:
    """
    stringtie --merge ${gtfs.flatten().join(' ')} -G $gtf -e -F -T -o merged.gtf
    """
  }

bam_stringtie2
    .combine(merged_gtf)
    .into {stringtie_input}

process stringtie2 {
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/stringtie2", mode: 'copy'

    input:
    set sampleid, file(bamfiles), file(mergedgtf) from stringtie_input

    output:
    file '*_transcripts.gtf' into final_gtf_list
    file '*.gene_abund.txt' into gene_abund
    file '*.cov_refs.gtf'
    stdout into stringtie_log

    script:
    """
    stringtie $bamfiles \\
    -o ${bamfiles}_transcripts.gtf \\
    -v \\
    -G $mergedgtf \\
    -A ${bamfiles}.gene_abund.txt \\
    -C ${bamfiles}.cov_refs.gtf \\
    -e \\
    -b ${bamfiles}_ballgown
    """
  }


  // Make gene count and transcript count matrix
  process aggregate_counts {
     publishDir "${params.outdir}/counts_files", mode: "copy"

     input:
     val abund_file_list from gene_abund.toList()

     output:
     file "fpkm.csv" into fpkm_counts
     file "tpm.csv" into tpm_counts

     script:
     String file_list = abund_file_list.flatten().join(' ')
     """
     python $PWD/aggrcounts.py $file_list
     """
  }

}

// MULTIQC
process multiqc {
  publishDir "${params.outdir}/multiqc_files", mode: 'copy'

  input:
  file ('fastqc/*')     from fastqc_results.flatten().toList()
  file ('trimgalore/*') from trimgalore_results.flatten().toList()
  file ('alignment/*')  from alignment_logs.flatten().toList()
  file ('rseqc/*')      from gene_coverage_results.flatten().toList()
  file ('rseqc/*')      from junction_annotation_results.flatten().toList()
  file ('stringtie/*')  from stringtie_log.flatten().toList()
  file ('counts/*')     from fpkm_counts.flatten().toList()

  output:
  file "*multiqc_report.html"
  file "*multiqc_data"

  script:
  println fastqc_results
  println fastqc_results.flatten().toList()
  """
  multiqc -f ${params.outdir}
  """
}

// ---------------------------------------------//

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}
