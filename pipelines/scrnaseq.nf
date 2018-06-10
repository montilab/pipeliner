#!/usr/bin/env nextflow

version = 2.6

// ---------------------------------------------//
//                 COPY CONFIG                  //
// ---------------------------------------------//

/*
Makes a copy of the configuration file used each run
and stores a timestamped copy in /history directory
*/

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
//             LOAD / VALIDATE FILES            //
// ---------------------------------------------//

/*
The pipeline can start from fastq read files or
start from alignment .bam files
*/

if (!params.bams.use_existing) {
  // Starting from .fastq files
  if (params.paired) {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_check; reads_preqc; reads_fastqc; reads_whitelist; reads_extract}
  } else {
    log.info "ERROR: Reads must be paired"
    System.exit(0)
  }
  /*
  ** Check reads
  */
  if (!params.skip.check_reads) {
    process check_reads {
      cache "deep"; tag "$sampleid"

      input:
      set sampleid, reads from reads_check

      script:
      template 'checks/fq/paired.sh'
    }
  }
} else {
  // Starting from .bam files
  Channel
    .fromPath(params.bams.path).splitCsv(header: true)
    .map {row -> [row.Sample_Name, file(row.Alignment)]}
    .ifEmpty {error "File ${params.alignments} not parsed properly"}
    .into {bam_files; bam_counts; bam_rseqc;}
}

/*
The pipeline requires the following under various conditions
- Annotation file is always required regardless of workflow
- Reference file is only required when building an index to align fastq files
*/

if (params.gtf) {
  Channel
    .fromPath(params.gtf).ifEmpty {exit 1, "GTF annotation file not found: ${params.gtf}"}
    .toList().into {gtf_indexing; gtf_mapping; gtf_parsing;}
}

// Either build index or load existing index
if (!params.bams.use_existing) {
  if (!params.index.use_existing) {
    // Build index if not using bams and no index exists
    if (params.fasta) {
      Channel
        .fromPath(params.fasta).ifEmpty {exit 1, "FASTA reference file not found: ${params.fasta}"}
        .toList().set {fasta_indexing;}
    }
    /*
    ** Build index with HISAT
    */
    if (params.aligner == 'hisat') {
      process hisat_indexing {
        cache "deep"; tag "$fasta"
        publishDir path: "${params.outdir}/${params.aligner}_index", mode: 'copy'

        input:
        file fasta from fasta_indexing
        file gtf from gtf_indexing

        output:
        file "index/*.ht2" into hisat_index

        script:
        template 'hisat/indexing.sh'
      }
    }
    /*
    ** Build index with STAR
    */
    if (params.aligner == 'star') {
      process star_indexing {
        cache "deep"; tag "$fasta"
        publishDir path: "${params.outdir}/${params.aligner}_index", mode: 'copy'

        input:
        file fasta from fasta_indexing
        file gtf from gtf_indexing

        output:
        file "index" into star_index

        script:
        template 'star/indexing.sh'
      }
    }
  } else {
    hisat_index = Channel.fromPath("${params.index.path}*")
    star_index = Channel.fromPath("${params.index.path}/*")
  }
}

// -----------------------------------------------------------------------------
//                BEGIN PIPELINE
// -----------------------------------------------------------------------------
if (!params.bams.use_existing) {
  /*
  ** Fastqc
  */
  if (!params.skip.fastqc) {
    process fastqc {
      cache "deep"; tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/fastqc", mode: 'copy'

      input:
      set sampleid, reads from reads_fastqc

      output:
      set sampleid, '*_fastqc.{zip,html}' into fastqc_results

      script:
      template 'fastqc/paired.sh'
    }
  } else {
    Channel.from(false).set {fastqc_results}
  }

  /*
  ** Whitelist
  */
  process whitelist {
    cache "deep"; tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/whitelisted", mode: 'copy'

    input:
    set sampleid, reads from reads_whitelist

    output:
    set sampleid, 'whitelist.txt' into whitelist

    script:
    template 'whitelist.sh'
  }
  
  /*
  ** Extract
  */
  process extract {
    cache "deep"; tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/extracted", mode: 'copy'

    input:
    set sampleid, reads from reads_extract
    file w from whitelist

    output:
    set sampleid, '*extracted.fastq.gz' into extracted_reads

    script:
    template 'extract.sh'
  }

  /*
  ** Mapping
  */
  if (params.aligner == "hisat") {
    process hisat_mapping {
      cache "deep"; tag "$sampleid"
      publishDir "${params.outdir}/${params.aligner}_bams", mode: "copy",
      saveAs: {filename ->
                if (filename.indexOf("log") > 0) "logs/$filename"
                else if (params.bams.save) filename
                else null}

      input:
      file index from hisat_index.collect()
      file gtf from gtf_mapping
      set sampleid, file (reads:"*") from extracted_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_rseqc
      set sampleid, '*.log' into mapping_results

      script:
      if (!params.index.use_existing) {
        index_base = index[0].toString() - ~/.\d.ht2/
      } else {
        index_base = params.index.path
      }

      template 'hisat/mapping/paired.merged.sh'
    }
  }
  if (params.aligner == "star") {
    process star_mapping {
      cache "deep"; tag "$sampleid"
      publishDir "${params.outdir}/${params.aligner}_bams", mode: "copy",
       saveAs: {filename ->
                if (filename.indexOf("Log") > 0) "logs/$filename"
                else if (params.bams.save) "${sampleid}.bam"
                else null}

      input:
      file index from star_index.collect()
      file gtf from gtf_mapping
      set sampleid, file (reads:"*") from extracted_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_rseqc
      set sampleid, '*.out' into mapping_results

      script:
      if (params.index.use_existing) {
        index = params.index.path
      }
      template "star/mapping.merged.sh"
    }
  }
} else {
  Channel.from(false).set {fastqc_results}
  Channel.from(false).set {mapping_results}
}

if (!params.skip.rseqc) {
  /*
  ** Bamqc
  */
  process gtftobed {
    cache "deep"; tag "$gtf"
    publishDir "${params.outdir}/rseqc", mode: 'copy'

    input:
    file gtf from gtf_mapping

    output:
    file "*.bed" into bed_files

    script:
    """
    python $PWD/scripts/rseqc/gfftobed.py ${gtf} > ${gtf}.bed
    """
  }
  process rseqc {
    cache "deep"
    tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/rseqc", mode: 'copy'

    input:
    set sampleid, file(bamfiles) from bam_rseqc
    file bed from bed_files

    output:
    file ("*.{txt,pdf,r,xls}")
    file ("*summary.txt") into rseqc_results
    
    script:
    template 'rseqc.sh'
  }
} else {
  Channel.from(false).set {rseqc_results}
}

if (!params.skip.counting) {
  if (params.quantifier == 'featurecounts') {
    /*
    ** Counting
    */
    process counting {
      cache "deep"; tag "$sampleid"
      publishDir "${params.outdir}/${sampleid}/${params.quantifier}", mode: 'copy'

      input:
      set sampleid, file(bamfiles) from bam_counts
      file gtf from file(params.gtf)

      output:
      file '*'
      file '*counts.txt' into counts
      stdout into quant_results

      script:
      template 'featurecounts.umitools.sh'
    }
  }
} else {
  Channel.from(false).set {quant_results}
}

if (!params.skip.multiqc) {
  /*
  ** Multiqc
  */
  process multiqc {
    cache "deep"
    publishDir "${params.outdir}/reports/post_pipeliner", mode: 'copy'

    input:
    file ('*') from fastqc_results.flatten().toList()
    file ('*') from mapping_results.flatten().toList()
    file ('*') from rseqc_results.flatten().toList()
    file ('*') from quant_results.flatten().toList()

    output:
    file "*report.html"
    file "*data"

    script:
    template 'multiqc/basic.sh'
  }
}

// -----------------------------------------------------------------------------

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}