#!/usr/bin/env nextflow

version = 2.3

// ---------------------------------------------//
//                VALIDATE CONFIG               //
// ---------------------------------------------//

cmd = workflow.commandLine
def range = 0..cmd.length()-6
config = ""
for (i in range) { 
  if (cmd[i..i+5] == ".nf -c") {
      config = cmd[i+7..cmd.length()-1]
      def src = new File("${config}")
      def dst = new File("${workflow.start}.config")
      dst << src.text
      break
  }
}
if (!config) {
  log.info "WARN: Copy of configuration file not made"
  log.info "Please specify a configuration file (e.g. -c nextflow.config)"
}

// ---------------------------------------------//
//                VALIDATE FILES                //
// ---------------------------------------------//

// Read and map reads with samples using csv file
if (!params.from_bam) {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_fastqc; reads_cutadapt}
} else {
  Channel
    .fromPath(params.alignments).splitCsv(header: true)
    .map {row -> [row.Sample_Name, file(row.Alignment)]}
    .ifEmpty {error "File ${params.alignments} not parsed properly"}
    .into {bam_files; bam_counts; bam_stringtie1; bam_stringtie2; bam_rseqc}
}

if (params.gtf) {
  Channel
    .fromPath(params.gtf).ifEmpty {exit 1, "GTF annotation file not found: ${params.gtf}"}
    .toList().into {gtf_indexing; gtf_mapping; gtf_stringtie;}
}

log.info " P I P E L I N E R  ~  v${version}"
log.info "===================================="
log.info "Reads          : ${params.reads}"
log.info "FASTA          : ${params.fasta}"
log.info "Annotation     : ${params.gtf}"
log.info "Input Dir      : ${params.indir}"
log.info "Output Dir     : ${params.outdir}"

log.info "===================================="
log.info "Aligner        : ${params.aligner}"
log.info "Save Reference : ${params.save_reference}"
log.info "Save Temporary : ${params.save_temps}"

log.info "===================================="
log.info "Current user  : $USER"
log.info "Current home  : $HOME"
log.info "Current path  : $PWD"
log.info "===================================="

// ---------------------------------------------//
//                BEGIN PIPELINE                //
// ---------------------------------------------//

if (!params.from_bam) {
  
  // FASTQC
  process fastqc {
    cache params.caching
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/fastqc", mode: 'copy'

    input:
    set sampleid, reads from reads_fastqc

    output:
    set sampleid, '*_fastqc.{zip,html}' into fastqc_results

    script:
    template 'fastqc.sh'
  }

  // CUTADAPT
  process cutadapt {
    cache params.caching
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/cutadapt", mode: 'copy'

    input:
    set sampleid, reads from reads_cutadapt

    output:
    set sampleid, '*.fq' into trimmed_reads
    set sampleid, '*.log' into cutadapt_log

    script:
    template 'cutadapt.sh'
  }

  // STAR - BUILD INDEX
  if (params.aligner == 'star' && file(params.fasta)) {
    process star_indexing {
      cache params.caching
      tag "$fasta"
      publishDir path: "${params.outdir}/star_files", saveAs: {params.save_reference ? it : null}, mode: 'copy'

      input:
      file fasta from file(params.fasta)
      file gtf from gtf_indexing

      output:
      file "star_index" into star_index

      script:
      template 'star_indexing.sh' 
    }
  }

  // STAR - MAP READS
  if (params.aligner == 'star') {
    process star_mapping {
      cache params.caching    
      tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/star", mode: 'copy'

      input:
      file index from star_index
      file gtf from gtf_mapping
      set sampleid, file (reads:'*') from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_rseqc, bam_htseq
      set sampleid, '*.out' into alignment_logs
      set sampleid, '*.tab' into alignment_tab

      script:
      template 'star_mapping.sh' 
    }
  }
}

// RSEQC
def num_bams
bam_counts.count().subscribe{num_bams=it}
process rseqc {
  cache params.caching
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/rseqc", mode: 'copy'

  input:
  set sampleid, file(bamfiles) from bam_rseqc

  output:
  file ('*.bam_stats') into bam_stats_results

  script:
  template 'rseqc.sh'
}

// HTSEQ
process htseq {
  cache params.caching    
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/htseq", mode: 'copy'

  input:
  set sampleid, file(bamfiles) from bam_htseq
  file gtf from file(params.gtf)

  output:
  file '*counts.txt' into counts
  stdout into quantification

  script:
  template 'htseq.sh'
}

// MULTIQC
process multiqc {
  cache params.caching
  publishDir "${params.outdir}/multiqc_files", mode: 'copy'

  input:
  file ('fastqc/*')   from fastqc_results.flatten().toList()
  file ('cutadapt/*') from cutadapt_log.flatten().toList()
  file ('rseqc/*')    from bam_stats_results.flatten().toList()
  file ('star/*')     from alignment_logs.flatten().toList()
  file ('htseq/*')    from quantification.flatten().toList()

  output:
  file "*multiqc_report.html"
  file "*multiqc_data"

  script:
  template 'multiqc.sh'
}

// ---------------------------------------------//

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}