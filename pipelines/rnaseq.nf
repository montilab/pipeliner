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
      .into {reads_check; reads_preqc; reads_trimming;}
  } else {
    Channel
      .fromPath(params.reads).splitCsv(header: true)
      .map {row -> [row.Sample_Name, [row.Read1]]}
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into {reads_check; reads_preqc; reads_trimming;}
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
      if (params.paired){
        template 'checks/fq/paired.sh'
      } else {
        template 'checks/fq/single.sh'
      }
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
//              PRE-PIPELINE QC
// -----------------------------------------------------------------------------
if (!params.bams.use_existing) {
  if (!params.skip.pre_qc) {
    process pre_fastqc {
      cache "deep"; tag "$sampleid"
      publishDir "${params.outdir}/pre_fastqc", mode: 'copy'

      input:
      set sampleid, reads from reads_preqc

      output:
      set sampleid, '*_fastqc.{zip,html}' into preqc_results

      script:
      if (params.paired){
        template 'fastqc/paired.sh'
      } else {
        template 'fastqc/single.sh'
      }
    }
    process pre_multiqc {
      cache "deep"
      publishDir "${params.outdir}/reports/pre_pipeliner", mode: 'copy'

      input:
      file ('*') from preqc_results.flatten().toList()

      output:
      file "*data"
      file "*report.html"

      script:
      """
      multiqc ${params.outdir}/pre_fastqc
      """
    }
  }
}

// -----------------------------------------------------------------------------
//                BEGIN PIPELINE
// -----------------------------------------------------------------------------
if (!params.bams.use_existing) {
  /*
  ** Trimming
  */
  process trim_galore {
    cache "deep"; tag "$sampleid"
    publishDir "${params.outdir}/${sampleid}/trimgalore", mode: 'copy',
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
      set sampleid, file (reads:"*") from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_rseqc
      set sampleid, '*.log' into mapping_results

      script:
      if (!params.index.use_existing) {
        index_base = index[0].toString() - ~/.\d.ht2/
      } else {
        index_base = params.index.path
      }

      if (params.paired){
        template "hisat/mapping/paired.sh"
      } else {
        template "hisat/mapping/single.sh"
      }
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
      set sampleid, file (reads:"*") from trimmed_reads

      output:
      set sampleid, file("*.bam") into bam_files, bam_counts, bam_rseqc
      set sampleid, '*.out' into mapping_results

      script:
      if (params.index.use_existing) {
        index = params.index.path
      }
      template "star/mapping.sh"
    }
  }
// -----------------------------------------------------------------------------
} else {
  Channel.from(false).set {trim_galore_results}
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
  if (params.quantifier == 'htseq' || params.quantifier == 'stringtie' || params.quantifier == 'featurecounts') {
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
      if (params.quantifier == 'htseq') {
        template 'htseq.sh'
      } else if (params.quantifier == 'stringtie') {
        template 'stringtie.sh'
      } else if (params.quantifier == 'featurecounts') {
        template 'featurecounts.sh'
      }
    }
  }
  /*
  ** Matrix
  */
  process expression_matrix {
    cache "deep"    
    publishDir "${params.outdir}/expression_matrix", mode: 'copy'

    input:
    val counts_files from counts.toList()

    output:
    file '*.txt'
    file '*expression_matrix.txt' into matrix_parsing, expression_matrix

    script:
    String files = counts_files.flatten().join(' ')
    if (!params.skip.eset) {
      """
      python $PWD/scripts/expression/create_matrix.py -p ${params.phenotypes} ${params.quantifier} ${files}
      """
    } else {
      """
      python $PWD/scripts/expression/create_matrix.py ${params.quantifier} ${files}
      """
    }
  }
} else {
  Channel.from(false).set {quant_results}
}

/*
** Eset
*/
if (!params.skip.eset) {
  process expression_features {
    cache "deep"    
    publishDir "${params.outdir}/expression_matrix", mode: 'copy'

    input:
    file gtf from gtf_parsing
    file matrix from matrix_parsing

    output:
    file 'features.txt' into features

    script:
    """
    python $PWD/scripts/expression/parse_gtf.py ${gtf} ${matrix}
    """     
  }
  process expression_set {
    cache "deep"    
    publishDir "${params.outdir}/expression_set", mode: 'copy'

    input:
    file matrix from expression_matrix
    file features from features
    output:
    file '*expression_set.rds' into expression_set

    script:
    """
    Rscript $PWD/scripts/expression/create_expressionset.R ${matrix} ${features} ${params.phenotypes}
    """
  }
}

if (!params.skip.multiqc) {
  /*
  ** Multiqc
  */
  process multiqc {
    cache "deep"
    publishDir "${params.outdir}/reports/post_pipeliner", mode: 'copy'

    input:
    file ('*') from trim_galore_results.flatten().toList()
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