#!/usr/bin/env nextflow

version = 2.6

if (params.paired) {
  Channel
    .fromPath(params.reads).splitCsv(header: true)
    .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
    .ifEmpty {error "File ${params.reads} not parsed properly"}
    .set {reads_trimming}
} else {
  Channel
    .fromPath(params.reads).splitCsv(header: true)
    .map {row -> [row.Sample_Name, [row.Read1]]}
    .ifEmpty {error "File ${params.reads} not parsed properly"}
    .set {reads_trimming}
}

process trim_galore {
  cache "deep"; tag "$sampleid"
  publishDir "${params.outdir}/trimgalore/${sampleid}", mode: 'copy',
  saveAs: {filename ->
            if (filename.indexOf("trimming_report.txt") > 0) filename
            else if (filename.indexOf("fastqc.html") > 0) filename
            else if (filename.indexOf("fastqc.zip") > 0) filename
            else if (params.fastqs.save) filename
            else null}

  input:
  set sampleid, reads from reads_trimming

  output:
  set sampleid, '*fq.gz' into trimmed_reads
  set sampleid, '*' into trim_galore_results

  script:
  if (params.trim_galore.custom_adaptors) {
    if (params.paired){
        template 'trim_galore/custom/paired.sh'
    } else {
      template 'trim_galore/custom/single.sh'
    }
  } else {
    if (params.paired){
        template 'trim_galore/default/paired.sh'
    } else {
      template 'trim_galore/default/single.sh'
    }        
  }
}

if (!params.skip.multiqc) {
  process multiqc {
    cache "deep"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file ('*') from trim_galore_results.flatten().toList()

    output:
    file "*report.html"
    file "*data"

    script:
    template 'multiqc/basic.sh'
  }
}

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
  def subject = 'Barebones Pipeline Execution'
  def recipient = (params.email)
    ['mail', '-s', subject, recipient].execute() << """
    Pipeline Summary
    ---------------------------
    Timestamp    : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit Status  : ${workflow.exitStatus}
    Error Report : ${workflow.errorReport ?: '-'}
    """
}