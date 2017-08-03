#!/usr/bin/env nextflow

// Pipeline version
version = 1.1

// Configurable variables
params.fasta = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data/genome.index.fa'
params.gtf   = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data/genome.bed.gff'
params.reads_file = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data/ggal_alpha.csv'
params.outdir = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/results/results'
params.bed = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data/ggal_refseq.bed'
params.aligner = 'star'
params.starindex = false
params.paired = true
params.saveReference = true

// Read and Map reads with samples using csv file
if (params.paired) {
    Channel
        .fromPath(params.reads_file)
        .splitCsv(header: true)
        .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
        .ifEmpty { error "File ${params.reads_file} not parsed properly" }
        .into {read_files_fastqc; read_files_trimming}
} else {
    Channel
        .fromPath(params.reads_file)
        .splitCsv(header: true)
        .map {row -> [row.Sample_Name, [row.Read]]}
        .ifEmpty { error "File ${params.reads_file} not parsed properly" }
        .into {read_files_fastqc; read_files_trimming}
}


// Validate inputs
if (params.starindex && params.aligner == 'star'){
    starindex = Channel
        .fromPath(params.starindex)
        .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
        .toList()
}

if (params.gtf){
    Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .toList()
    .into { gtf_makeSTARindex; gtf_bowtieindex;gtf_star;gtf_bowtie;gtf_stringtieFPKM;  }
}

if( params.bed ){
    bed = Channel
        .fromPath(params.bed)
        .ifEmpty { exit 1, "BED annotation file not found: ${params.bed}" }
        .toList()

}
fasta = file(params.fasta)
gtf   = file(params.gtf)
bed = file(params.bed)


log.info " P I P E L I N E R  ~  v${version}"
log.info "===================================="
log.info "Reads         : ${params.reads}"
log.info "Genome        : ${params.genome}"
log.info "FASTA         : ${params.fasta}"
log.info "Annotation    : ${params.gtf}"
log.info "Output dir    : ${params.outdir}"

        if(params.aligner == 'star')
        {
        log.info "Aligner        : STAR"
                if (params.starindex)
                        log.info "STAR Index: ${params.starindex}"

        }

log.info "Current home  : $HOME"
log.info "Current user  : $USER"
log.info "Current path  : $PWD"
log.info "===================================="

logParams(params, "nextflow_paramters.txt")

def logParams(p, n) {
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
  set sampleid, reads from read_files_fastqc

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

// TRIM GALORE
process trim_galore {
  tag "$sampleid"
  cache 'deep'
  publishDir "${params.outdir}/$sampleid/trim_galore", mode: 'copy'

  input:
  set sampleid, reads from read_files_trimming

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

// STAR - BUILD INDEX
if (params.aligner == 'star' && !params.starindex && fasta) {
  process makeSTARindex {

    publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'
    tag "$fasta"
    cache 'deep'

    input:
    file fasta from fasta
    file gtf from gtf_makeSTARindex

    output:
    file "star" into starindex

    script:
    """
    mkdir star
    STAR \\
    --runMode genomeGenerate \\
    --runThreadN ${task.cpus} \\
    --sjdbGTFfile $gtf \\
    --sjdbOverhang 149 \\
    --genomeDir star/ \\
    --genomeFastaFiles $fasta
    """
  }
}

// STAR - MAP READS
if (params.aligner == 'star') {
  process star {
    tag "$sampleid"
    cache 'deep'
    publishDir "${params.outdir}/$sampleid/STAR", mode: 'copy'

    input:
    file index from starindex.first()
    file gtf from gtf_star.first()
    set sampleid, file (reads:'*') from trimmed_reads

    output:
    set sampleid, file("*.bam") into  bam_files, bam_count, bam_stringtieFPKM1, bam_stringtieFPKM2, bam_rseqc
    set sampleid, '*.out' into alignment_logs
    set sampleid, '*SJ.out.tab' into alignment_tab

    script:
    """
    prefix='${sampleid}.'
    STAR --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --readFilesIn $reads  \\
    --runThreadN ${task.cpus} \\
    --twopassMode Basic \\
    --outWigType bedGraph \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat \\
    --outFileNamePrefix \$prefix
    """
  }
}

// MULTIQC
process multiqc {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  input:
  file ('fastqc/*') from fastqc_results.flatten().toList()
  file ('trimgalore/*') from trimgalore_results.flatten().toList()
  file ('alignment/*') from alignment_logs.flatten().toList()
  //file ('stringtie/*') from stringtie_log.flatten().toList()
  //file('counts/*') from fpkm_counts.flatten().toList()
  //file ('rseqc/*') from gene_coverage_results.flatten().toList()
  //file ('rseqc/*') from junction_annotation_results.flatten().toList()

  output:
  file "*multiqc_report.html"
  file "*multiqc_data"

  script:
  """
  multiqc -f ${params.outdir}
  """
}


// ---------------------------------------------//

workflow.onComplete {
  println ( workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong." )
}