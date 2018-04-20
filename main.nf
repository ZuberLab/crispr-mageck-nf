#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
     crispr-mageck-nf
    ================================================================

     Statistical analysis of multiplexed CRISPR-Cas9 / shRNA screens
     using MAGeCK.

     Usage:
     nextflow run zuberlab/crispr-mageck-nf

     Options:
        --contrasts     Tab-delimited text file specifying the contrasts
                        to be analyzed. (Defaults to 'contrasts.txt')
                        The following columns are required:
                            - name: name of contrasts
                            - control: control samples (comma separated)
                            - treatment: treatment samples (comma separated)
                            - norm_method: normalization method
                            - fdr_method: multiple testing adjustment method
                            - lfc_method: method to combine guides / hairpins
                            - cnv_correction

        --counts        Tab-delimited test file containing the raw counts.
                        (Defaults to 'counts_mageck.txt')
                        This file must conform to the input requirements of
                        MAGeCK 0.5.7 (http://mageck.sourceforge.net)

        --outputDir     Directory name to save results to. (Defaults to
                        '03_stats')

     Profiles:
        standard        local execution
        singularity     local execution with singularity
        ii2             IMPIMBA2 cluster execution with singularity

     Docker:
     zuberlab/crispr-nf:latest

     Author:
     Jesse J. Lipp (jesse.lipp@imp.ac.at)

    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

Channel
    .fromPath( params.contrasts )
    .splitCsv(sep: '\t', header: true)
    .set { contrastsMageck }

Channel
    .fromPath( params.counts )
    .set { countsMageck }

Channel
    .fromPath( params.cnv )
    .set { cnvMageck }

process mageck {

    tag { parameters.name }

    publishDir path: "${params.outputDir}/${parameters.name}",
               mode: 'copy',
               overwrite: 'true',
               saveAs: {filename ->
                   if (filename.indexOf(".log") > 0) "$filename"
                   else if (filename.indexOf(".normalized.txt") > 0) "$filename"
                   else null
               }
    input:
    val(parameters) from contrastsMageck
    each file(counts) from countsMageck
    each file(cnv) from cnvMageck

    output:
    set val("${parameters.name}"), file('*.sgrna_summary.txt'), file('*.gene_summary.txt') into resultsMageck
    file('*.log') into logsMageck
    file('*.normalized.txt') into normalizedMageck

    script:
    rra_params = params.min_rra_window > 0 ? "--additional-rra-parameters '-p ${params.min_rra_window}'" : ''
    cnv_file = cnv.exists() & parameters.cnv_correction != '' ? "--cnv-norm ${cnv}" : ""
    cnv_cellline = cnv.exists() & parameters.cnv_correction != '' ? "--cell-line ${parameters.cnv_correction}" : ""
    """
    prefilter_counts.R \
        ${counts} \
        ${parameters.control} \
        ${params.min_count} > counts_filtered.txt

    mageck test \
        --output-prefix ${parameters.name} \
        --count-table counts_filtered.txt \
        --control-id ${parameters.control} \
        --treatment-id ${parameters.treatment} \
        --norm-method ${parameters.norm_method} \
        --adjust-method ${parameters.fdr_method} \
        --gene-lfc-method ${parameters.lfc_method} \
        --normcounts-to-file \
        ${rra_params} \
        ${cnv_file} \
        ${cnv_cellline}
    """
}

process postprocess {

    tag { name }

    publishDir path: "${params.outputDir}/${name}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(name), file(guides), file(genes) from resultsMageck

    output:
    set val(name), file('*_stats.txt') into processedMageck
    file('*.pdf') into qcMageck

    script:
    """
    postprocess_mageck.R ${guides} ${genes}
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
