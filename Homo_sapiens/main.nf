#!/usr/bin/env nextflow

// Pipeline version
version = 2.0

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

log.info " P I P E L I N E R  ~  v${version}"
log.info "===================================="
log.info "Reads         : ${params.reads}"
log.info "FASTA         : ${params.fasta}"
log.info "Annotation    : ${params.gtf}"
log.info "Input Dir     : ${params.indir}"
log.info "Output Dir    : ${params.outdir}"
log.info "Aligner       : ${params.aligner}"
if (params.index) {
  log.info "Index       : ${params.index}"
}
log.info "Current user : $USER"
log.info "Current home : $HOME"
log.info "Current path : $PWD"
log.info "===================================="

// ---------------------------------------------//
//                BEGIN PIPELINE                //
// ---------------------------------------------//

// FASTQC
process fastqc {
  tag "$sampleid"
  publishDir "${params.outdir}/$sampleid/fastqc", mode: 'copy'
  if (params.enable_resume) {
    cache 'deep'
  }

  input:
  set sampleid, reads from reads_fastqc

  output:
  set sampleid, '*_fastqc.{zip,html}' into fastqc_results

  script:
  """
  fastqc --quiet --outdir . ${reads[0]} ${reads[1]}
  """
}

// ---------------------------------------------//

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}
