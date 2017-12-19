#!/usr/bin/env nextflow

// Pipeline version
version = 2.4

// ---------------------------------------------//
//                  COPY CONFIG                 //
// ---------------------------------------------//

cmd = workflow.commandLine
def range = 0..cmd.length()-6
config = ""
for (i in range) { 
  if (cmd[i..i+5] == ".nf -c") {
      config = cmd[i+7..cmd.length()-1]
      def folder = new File('history')
      if( !folder.exists() ) {
        folder.mkdirs()
      }
      def filename = "${workflow.start}".replaceAll("[ :]", "_");
      def src = new File("${config}")
      def dst = new File("history/${filename}.config")
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

if (!params.from_bam) {
  if (params.paired) {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_check; reads_fastqc; reads_trimgalore;}
  } else {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read1]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_check; reads_fastqc; reads_trimgalore;}
  }
  if (!params.check_reads.skip) {
    process check_reads {
      cache params.caching
      tag "$sampleid"

      input:
      set sampleid, reads from reads_check

      script:
      if (params.paired){
        """
        python $PWD/scripts/quality/check_reads.py ${reads[0]} ${reads[1]}
        """
      } else {
        """
        python $PWD/scripts/quality/check_reads.py ${reads[0]}
        """
      }
    }
  }
} else {
  Channel
    .fromPath(params.alignments).splitCsv(header: true)
    .map {row -> [row.Sample_Name, file(row.Alignment)]}
    .ifEmpty {error "File ${params.alignments} not parsed properly"}
    .into {bam_files; bam_counts; bam_stringtie1; bam_stringtie2; bam_rseqc;}
}

if (params.fasta) {
  Channel
    .fromPath(params.fasta).ifEmpty {exit 1, "FASTA reference file not found: ${params.fasta}"}
    .toList().into {fasta_indexing;}
}

if (params.gtf) {
  Channel
    .fromPath(params.gtf).ifEmpty {exit 1, "GTF annotation file not found: ${params.gtf}"}
    .toList().into {gtf_indexing; gtf_mapping;}
}

log.info " P I P E L I N E R  ~  v${version}"

log.info "===================================="
log.info "Reads          : ${params.reads}"
log.info "Reference      : ${params.fasta}"
log.info "Annotation     : ${params.gtf}"
log.info "Input Dir      : ${params.indir}"
log.info "Output Dir     : ${params.outdir}"

log.info "===================================="
if (params.paired) {
  log.info "Read Type      : paired-end"
} else {
  log.info "Read Type      : single-end"
}
log.info "Aligner        : ${params.aligner}"
if (params.pre_indexed) {
  log.info "Index        : ${params.index}"
}
log.info "Quantifier     : ${params.quantifier}"
log.info "Save Reference : ${params.save_reference}"
log.info "Save Temporary : ${params.save_temps}"

log.info "===================================="
log.info "Current user  : $USER"
log.info "Current home  : $HOME"
log.info "Current path  : $PWD"
log.info "===================================="


// -----------------------------------------------------------------------------
//              PRE-PIPELINE FASTQC
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
//                BEGIN PIPELINE
// -----------------------------------------------------------------------------

if (!params.from_bam) {

// -----------------------------------------------------------------------------
// TRIMGALORE - Trims low quality or contaminated reads
// -----------------------------------------------------------------------------
  if (!params.trim_galore.skip) {
    process trim_galore {
      cache params.caching
      tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/trimgalore", mode: 'copy'

      input:
      set sampleid, reads from reads_trimgalore

      output:
      set sampleid, '*fq.gz' into trimmed_reads, trimmed_reads_fastqc
      set sampleid, '*.txt' into trimgalore_results

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
  } else {
    Channel.from().set {trimgalore_results}
  }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// FASTQC - Quality control on read files
// -----------------------------------------------------------------------------
  if (!params.fastqc.skip) {
    process fastqc {
      cache params.caching
      tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/fastqc", mode: 'copy'

      input:
      set sampleid, file (reads:'*') from trimmed_reads_fastqc

      output:
      set sampleid, '*_fastqc.{zip,html}' into fastqc_results

      script:
      if (params.paired){
          template 'fastqc/paired.sh'
      } else {
        template 'fastqc/single.sh'
      }
    }
  } else {
    Channel.from().set {fastqc_results}
  }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// HISAT - A method for alignment reads to a reference and annotation file
// -----------------------------------------------------------------------------
  if (params.aligner == 'hisat' && !params.index && file(params.fasta)) {
    process hisat_indexing {
      cache params.caching
      tag "$fasta"
      publishDir path: "${params.outdir}/hisat_files", saveAs: {params.save_reference ? it : null}, mode: 'copy'

      input:
      file fasta from fasta_indexing
      file gtf from gtf_indexing

      output:
      file hisat_index into hisat_index

      script:
      template 'hisat/indexing.sh'
    }
  }
  if (params.aligner == 'hisat') {
    process hisat_mapping {
      cache params.caching    
      tag "$sampleid"
      scratch = true

      input:
      file index from hisat_index
      file gtf from gtf_mapping
      set sampleid, file (reads:'*') from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_stringtie1, bam_stringtie2, bam_rseqc
      set sampleid, '*.log' into alignment_logs

      script:
      if (params.paired){
          template 'hisat/mapping/paired.sh'
      } else {
        template 'hisat/mapping/single.sh'
      }
    }
  }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// STAR - A method for alignment reads to a reference and annotation file
// -----------------------------------------------------------------------------
  if (params.aligner == 'star' && !params.index && file(params.fasta)) {
    process star_indexing {
      cache params.caching
      tag "$fasta"
      publishDir path: "${params.outdir}/star_files", saveAs: {params.save_reference ? it : null}, mode: 'copy'

      input:
      file fasta from fasta_indexing
      file gtf from gtf_indexing

      output:
      file "star_index" into star_index

      script:
      template 'star/indexing.sh'
    }
  }
  if (params.aligner == 'star') {
    process star_mapping {
      cache params.caching    
      tag "$sampleid"
      scratch = true

      input:
      file index from star_index
      file gtf from gtf_mapping
      set sampleid, file (reads:'*') from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_stringtie1, bam_stringtie2, bam_rseqc
      set sampleid, '*.out' into alignment_logs

      script:
      template 'star/mapping.sh'
    }
  }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Writes alignment logs to a specific logging directory
// -----------------------------------------------------------------------------
  process write_alignment_logs {
    cache params.caching    
    tag "$sampleid"
    publishDir "${params.outdir}/logs/${params.aligner}", mode: 'copy' 

    input:
    set sampleid, file (logs:'*') from alignment_logs

    output:
    set sampleid, '*' into alignment_results

    script:
    if (params.aligner == 'hisat') {
      """
      cp ${logs} ${sampleid}.log
      """
    } else if (params.aligner == 'star') {
      """
      cp ${logs[0]} ${sampleid}.Log.final.out
      """
    }
  }
// -----------------------------------------------------------------------------

 
// -----------------------------------------------------------------------------
// Saves alignment files and writes paths to a csv file
// -----------------------------------------------------------------------------
  if (params.save_alignments){
    process save_alignments {
      cache params.caching    
      tag "$sampleid"
      publishDir "${params.outdir}/alignments", mode: 'copy'

      input:
      set sampleid, file (alignments:'*') from bam_files

      output:
      set sampleid, file("*.bam") into saved_alignments
      stdout into saving_out

      script:
      """
      cp ${alignments} ${sampleid}.bam
      """    
    }
    process clear_alignment_paths {
      input:
      file ('*') from saving_out.flatten().toList()

      output:
      stdout into clearing_out

      script:
      String outpath = "${params.outdir}/alignments/alignments.csv"
      """
      rm -rf ${params.outdir}/alignments/alignments.csv
      touch ${params.outdir}/alignments/alignments.csv
      echo Sample_Name,Alignment >> ${outpath}
      """   
    }
    process write_alignment_paths {
      cache params.caching    
      tag "$sampleid"
      publishDir "${params.outdir}/alignments", mode: 'copy'

      input:
      set sampleid, file (alignments:'*') from saved_alignments
      file ('*') from clearing_out.flatten().toList()

      script:
      String outpath = "${params.outdir}/alignments/alignments.csv"
      """
      echo ${sampleid},${params.outdir}/alignments/${sampleid}.bam >> ${outpath}
      """
    }
  }
// -----------------------------------------------------------------------------
} else {
  Channel.from().set {fastqc_results}
  Channel.from().set {trimgalore_results}
  Channel.from().set {alignment_results}
}


// -----------------------------------------------------------------------------
// RSEQC - Quality control on alignment files
// -----------------------------------------------------------------------------
if (!params.rseqc.skip) {
  process gtftobed {
    cache params.caching
    tag "$gtf"
    publishDir "${params.outdir}/bed_files", mode: 'copy'

    input:
    file gtf from gtf_mapping

    output:
    file "*.bed" into bed_files

    script:
    """
    python $PWD/scripts/rseqc/gfftobed.py $gtf > out.bed
    """
  }
  process rseqc {
    cache params.caching
    tag "$sampleid"
    publishDir "${params.outdir}/$sampleid/rseqc", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_rseqc
    file bed from bed_files

    output:
    file ("*.{txt,pdf,r,xls}") into rseqc_stats
    file ("*summary.txt") into rseqc_results
    
    script:
    template 'rseqc.sh'
  }
} else {
  Channel.from().set {rseqc_results}
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// HTSEQ - Quantification method for generation raw counts
// -----------------------------------------------------------------------------
if (params.quantifier == 'htseq') {
  process htseq {
    cache params.caching    
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/htseq", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_stringtie1
    file gtf from file(params.gtf)

    output:
    file '*counts.txt' into counts
    stdout into quant_results

    script:
    template 'htseq.sh'
  }
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// FEATURECOUNTS - Quantification method for generation raw counts
// -----------------------------------------------------------------------------
else if (params.quantifier == 'featurecounts') {
  process featurecounts {
    cache params.caching    
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/featurecounts", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_stringtie1
    file gtf from file(params.gtf)

    output:
    file '*counts.txt' into counts
    stdout into quant_results

    script:
    template 'featurecounts.sh'
  }
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// STRINGTIE - Quantification method for generation normalized counts
// -----------------------------------------------------------------------------
else {
  process stringtie {
    cache params.caching    
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/stringtie1", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_stringtie1
    file gtf from file(params.gtf)

    output:
    file '*.gtf' into gtf_list

    script:
    template 'stringtie1.sh'
  }

  process stringtie_merge {
    cache params.caching
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

  process stringtie_eb {
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/stringtie2", mode: 'copy'

    input:
    set sampleid, file(bamfiles), file(mergedgtf) from stringtie_input

    output:
    file '*.gtf'
    file '*.gene_abund.txt' into gene_abund
    file '*.txt' into quant_results

    script:
    template 'stringtie2.sh'
  }

  // AGGREGATE COUNTS
  process aggregate_counts {
    cache params.caching
    publishDir "${params.outdir}/counts", mode: "copy"

    input:
    val abundances from gene_abund.toList()

    output:
    file "fpkm.csv" into fpkm_results
    file "tpm.csv" into tpm_results

    script:
    String files = abundances.flatten().join(' ')
    """
    python $PWD/scripts/expression/normalize_counts.py $files
    """
  }
}

if (params.quantifier == 'htseq' || params.quantifier == 'featurecounts') {
  process expression_matrix {
    cache params.caching    
    publishDir "${params.outdir}/expression_matrix", mode: 'copy'

    input:
    val counts_files from counts.toList()

    output:
    file '*expression_matrix.txt' into expression_matrix

    script:
    String files = counts_files.flatten().join(' ')
    if (params.expression_set) {
      """
      python $PWD/scripts/expression/create_matrix.py -p ${params.phenotypes} ${params.quantifier} ${files}
      """
    } else {
      """
      python $PWD/scripts/expression/create_matrix.py ${params.quantifier} ${files}
      """
    }
  }
  if (params.expression_set) {
    process expression_set {
      cache params.caching    
      publishDir "${params.outdir}/expression_set", mode: 'copy'

      input:
      file matrix from expression_matrix

      output:
      file '*expression_set.rds' into expression_set

      script:
      """
      Rscript $PWD/scripts/expression/create_expressionset.R ${matrix} ${params.phenotypes} ${params.metadata}
      """
    }
  }
}

if (!params.multiqc.skip) {
  // MULTIQC
  process multiqc {
    cache params.caching
    publishDir "${params.outdir}/multiqc_files", mode: 'copy'

    input:
    file ('*') from trimgalore_results.flatten().toList()
    file ('*') from fastqc_results.flatten().toList()
    file ('*') from quant_results.flatten().toList()
    file ('*') from rseqc_results.flatten().toList()
    file ('*') from '${params.outdir}/logs/${params.aligner}'

    output:
    file "*report.html"
    file "*data"

    script:
    template 'multiqc.sh'
  }
}

// ---------------------------------------------//

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}