#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
        --contrasts      Tab-delimited text file specifying the contrasts
                         to be analyzed. (default: 'contrasts.txt')
                         The following columns are required:
                            - name: name of contrasts
                            - control: control samples (comma separated)
                            - treatment: treatment samples (comma separated)
                            - norm_method: normalization method
                            - fdr_method: multiple testing adjustment method
                            - lfc_method: method to combine guides / hairpins
                            - cnv_correction: cellline name
                            - filter : column to do count filtering on (default if empty: control)

        --counts         Tab-delimited text file containing the raw counts.
                         (default: 'counts_mageck.txt')
                         This file must conform to the input requirements of
                         MAGeCK 0.5.9 (http://mageck.sourceforge.net)

        --cnv            Tab-delimited text file containing the estimated
                         copy number per gene. Must be the same as in
                         'cnv_correction' in 'contrast' file. (default: 'cnv.txt')

        --min_count      Integer specifying the minimal control count to
                         consider the guide for further analysis. Guides with
                         lower counts in any control sample will be removed
                         prior to normalization for the given contrast.
                         (default: 50)

        --estimate_min_count_from_samples:  boolean indicating if min required counts for filtering
                                            should be estimated from the samples (> 5% of the mean count)

        --min_rra_window Float between 0 and 1 specifying the mininal fraction
                         of guides to consider for the RRA rank algorithm.
                         This overrides the dynamically determined threshold
                         and will use the same value for both negative and
                         positive selection. A value of 0 uses the dynamic
                         thresholding. (default: 0)

        --outputDir      Directory name to save results to. (Defaults to
                         '03_stats')

        --control_sgRNAs: txt file that contains sgRNA-ids to be used by mageck for normalization
                          and for generating the null distribution of RRA

     Profiles:
        standard         local execution
        apptainer        local execution with apptainer
        cbe              CBE cluster execution with apptainer

     Docker:
     zuberlab/crispr-nf:1.0

     Author:
     Florian Andersch (florian.andersch@imp.ac.at)

    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

// Import modules
include { MAGECK } from './modules/mageck'
include { POSTPROCESS } from './modules/postprocess'

contrastsMageck = Channel
    .fromPath(params.contrasts)
    .splitCsv(sep: '\t', header: true)

countsMageck = Channel.fromPath(params.counts)

cnvMageck = Channel.fromPath(params.cnv)

ctrlsgRNAs = Channel.fromPath(params.control_sgRNAs)

workflow {

    // Combine input channels
    ch_mageck = contrastsMageck
        .combine(countsMageck)
        .combine(cnvMageck)
        .combine(ctrlsgRNAs)

    // Run MAGeCK
    MAGECK(ch_mageck)
    
    ch_post_process = MAGECK.out.results
        .combine(ctrlsgRNAs)

    //  Postprocess results
    POSTPROCESS(ch_post_process)
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
