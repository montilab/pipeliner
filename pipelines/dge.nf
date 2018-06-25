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

Channel
  .fromPath(params.bams).splitCsv(header: true)
  .map {row -> [file(row.Alignment)]}
  .ifEmpty {error "File ${params.bams} not parsed properly"}
  .set {bam_counts;} 

/*
The pipeline requires the following under various conditions
- Annotation file is always required regardless of workflow
*/

if (params.gtf) {
  Channel
    .fromPath(params.gtf).ifEmpty {exit 1, "GTF annotation file not found: ${params.gtf}"}
    .toList().set {gtf_parsing;}
}

// -----------------------------------------------------------------------------
//                BEGIN PIPELINE
// -----------------------------------------------------------------------------

/*
** Counting
*/
process counting {
  publishDir "${params.outdir}/${params.quantifier}", mode: 'copy',
  saveAs: {filename ->
            if (filename.indexOf("txt.summary") > 0) filename
            else if (params.save.raw) filename
            else null}

  input:
  file(bamfiles) from bam_counts.collect()
  file gtf from file(params.gtf)

  output:
  file 'counts.raw.txt' into raw_counts
  file '*.summary' into quant_results

  script:
  feature_counts_sargs = ""
  if (params.paired) {
    feature_counts_sargs = feature_counts_sargs.concat("-p ")
  }
  if (params.remove_dups) {
    feature_counts_sargs = feature_counts_sargs.concat("--ignoreDup ")
  }
  template "feature_counts/vanilla.sh"
}

/*
** Matrix
*/
process expression_matrix {
  publishDir "${params.outdir}/expression_matrix", mode: 'copy'

  input:
  file counts from raw_counts

  output:
  file 'counts.matrix.txt' into matrix_counts

  script:
  """
  cut -f1,7- ${counts} | sed 1d > counts.matrix.txt
  """
}

/*
** Renaming
*/
process rename_samples {
  publishDir "${params.outdir}/expression_matrix", mode: 'copy'

  input:
  file matrix from matrix_counts

  output:
  file 'counts.matrix.renamed.txt'

  script:
  """
  python $PWD/scripts/expression/rename_samples.py ${matrix} ${params.bams}
  """
}

if (!params.skip.multiqc) {
  /*
  ** Multiqc
  */
  process multiqc {
    publishDir "${params.outdir}/reports/post_pipeliner", mode: 'copy'

    input:
    file ('*') from quant_results.flatten().toList()

    output:
    file "*report.html"
    file "*data"

    script:
    template "multiqc/vanilla.sh"
  }
}

// -----------------------------------------------------------------------------

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}