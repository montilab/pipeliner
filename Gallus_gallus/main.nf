#!/usr/bin/env nextflow

// Pipeline version
version = 2.1

// ---------------------------------------------//
//                VALIDATE FILES                //
// ---------------------------------------------//

// Read and map reads with samples using csv file
if (!params.from_bam) {
  if (params.paired) {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_fastqc; reads_trimgalore}
  } else {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_fastqc; reads_trimgalore}
  }
} else {
  Channel
    .fromPath(params.alignments).splitCsv(header: true)
    .map {row -> [row.Sample_Name, file(row.Alignment)]}
    .ifEmpty {error "File ${params.alignments} not parsed properly"}
    .into {bam_files; bam_counts; bam_stringtie1; bam_stringtie2}
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
log.info "Aligner        : ${params.aligner}"
if (params.pre_indexed) {
  log.info "Index        : ${params.enable_resume}"
}
log.info "===================================="
log.info "Resume Enabled : ${params.enable_resume}"
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
    cache params.enable_resume
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/fastqc", mode: 'copy'

    input:
    set sampleid, reads from reads_fastqc

    output:
    set sampleid, '*_fastqc.{zip,html}' into fastqc_results

    script:
    template 'fastqc.sh'
  }

  // TRIM GALORE
  process trim_galore {
    cache params.enable_resume
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/trimgalore", mode: 'copy'

    input:
    set sampleid, reads from reads_trimgalore

    output:
    set sampleid, '*fq.gz' into trimmed_reads
    set sampleid, '*trimming_report.txt' into trimgalore_results

    script:
    template 'trim_galore.sh'
  }

  // STAR - BUILD INDEX
  if (params.aligner == 'star' && !params.index && file(params.fasta)) {
    process star_indexing {
      cache params.enable_resume
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
      cache params.enable_resume    
      tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/star", mode: 'copy'

      input:
      file index from star_index
      file gtf from gtf_mapping
      set sampleid, file (reads:'*') from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_stringtie1, bam_stringtie2
      
      set sampleid, '*.out' into alignment_logs
      set sampleid, '*.tab' into alignment_tab

      script:
      template 'star_mapping.sh' 
    }
  }
}

// STRINGTIE
process stringetie {
  cache params.enable_resume    
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/stringtie1", mode: 'copy'

  input:
  set sampleid, file(bamfiles) from bam_stringtie1
  file gtf from file(params.gtf)

  output:
  file '*_transcripts.gtf' into gtf_list

  script:
  template 'stringtie1.sh'
}

process stringetie_merge {
  cache params.enable_resume
  publishDir "${params.outdir}/stringtiemerge", mode: 'copy'

  input:
  val gtfs from gtf_list.toList()
  file gtf from file(params.gtf)

  output:
  file 'merged.gtf' into merged_gtf

  script:
  template 'gtf_merging.sh'
}

bam_stringtie2
  .combine(merged_gtf)
  .set {stringtie_input}

process stringetie_eb {
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
  template 'stringtie2.sh'
}

// AGGREGATE COUNTS
process aggregate_counts {
  cache params.enable_resume
  publishDir "${params.outdir}/counts", mode: "copy"

  input:
  val abundances from gene_abund.toList()

  output:
  file "fpkm.csv" into fpkm_results
  file "tpm.csv" into tpm_results

  script:
  String files = abundances.flatten().join(' ')
  """
  python $PWD/scripts/aggregator.py $files
  """
}

if (!params.from_bam) {
  // MULTIQC
  process multiqc {
    cache params.enable_resume
    publishDir "${params.outdir}/multiqc_files", mode: 'copy'

    input:
    file ('fastqc/*')     from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('star/*')       from alignment_logs.flatten().toList()
    file ('stringtie/*')  from stringtie_log.flatten().toList()
    file ('counts/*')     from fpkm_results.flatten().toList()

    output:
    file "*multiqc_report.html"
    file "*multiqc_data"

    script:
    template 'multiqc.sh'
  }
}

// ---------------------------------------------//

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
  logger(params, "${params.outdir}/") // Requires outdir exists
}

def logger(p, dir) {
  filename = "${workflow.start}.log".replace(':','_')
  File file = new File(dir+filename)
  file.write "Work Directory: ${workflow.workDir}\n"
  file << "Command: ${workflow.commandLine}\n"
  file << "Start: ${workflow.start}\n"
  file << "End: ${workflow.complete}\n"
  file << "Duration: ${workflow.duration}\n"
  file << "Success: ${workflow.success}\n"
  file << "\nParameters\n"
  for (s in p) {
      file << "${s.key}: ${s.value}\n"
  } 
}
