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

        --legacy         Run Legacy MAGeCK 0.5.5

        --counts         Tab-delimited text file containing the raw counts.
                         (default: 'counts_mageck.txt')
                         This file must conform to the input requirements of
                         MAGeCK 0.5.7 (http://mageck.sourceforge.net)

        --cnv            Tab-delimited text file containing the estimated
                         copy number per gene. Must be the same as in
                         'cnv_correction' in 'contrast' file. (default: 'cnv.txt')

        --min_count      Integer specifying the minimal control count to
                         consider the guide for further analysis. Guides with
                         lower counts in any control sample will be removed
                         prior to normalization for the given contrast.
                         (default: 50)

        --min_rra_window Float between 0 and 1 specifying the mininal fraction
                         of guides to consider for the RRA rank algorithm.
                         This overrides the dynamically determined threshold
                         and will use the same value for both negative and
                         positive selection. A value of 0 uses the dynamic
                         thresholding. (default: 0)

        --outputDir      Directory name to save results to. (Defaults to
                         '03_stats')

     Profiles:
        standard         local execution
        singularity      local execution with singularity
        ii2              IMPIMBA2 cluster execution with singularity

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
    variance_params = parameters.variance_estimation ? "--variance-estimation-samples ${parameters.variance_estimation}" : ''
    rra_params = params.min_rra_window > 0 ? "--additional-rra-parameters '-p ${params.min_rra_window}'" : ''
    cnv_file = file(params.cnv).exists() & parameters.cnv_correction != '' ? "--cnv-norm ${cnv}" : ""
    cnv_cellline = file(params.cnv).exists() & parameters.cnv_correction != '' ? "--cell-line ${parameters.cnv_correction}" : ""

    control = parameters.filter == "" ? parameters.control : parameters.filter

    if( parameters.norm_method == "quantile" )

	    """
      IFS=',' read -ra CTRL <<< "${parameters.control}"
      ctrl_string=""
      for i in "\${CTRL[@]}"; do
        if [ "\${ctrl_string}" == "" ]; then
          ctrl_string_buff="\${i}_ref_${parameters.name}"
        else
          ctrl_string_buff=",\${i}_ref_${parameters.name}"
        fi
        ctrl_string="\${ctrl_string}\${ctrl_string_buff}"
      done

	    prefilter_counts.R \
	        ${counts} \
	        \${ctrl_string} \
          ${parameters.treatment} \
	        0 > counts_filtered.txt

	    quantile_normalize_counts.R \
	        counts_filtered.txt > counts_quantile_normalized.txt

	    VERSION=\$(mageck -v 2>&1 >/dev/null)

	    if [ \$VERSION = "0.5.5" ]; then

		    mageck test \
		        --output-prefix ${parameters.name} \
		        --count-table counts_quantile_normalized.txt \
		        --control-id \${ctrl_string_buff} \
		        --treatment-id ${parameters.treatment} \
		        --norm-method none \
		        --adjust-method ${parameters.fdr_method} \
		        --gene-lfc-method ${parameters.lfc_method} \
		        --normcounts-to-file

	    else
		    mageck test \
			        --output-prefix ${parameters.name} \
			        --count-table counts_quantile_normalized.txt \
			        --control-id \${ctrl_string_buff} \
			        --treatment-id ${parameters.treatment} \
			        --norm-method none \
			        --adjust-method ${parameters.fdr_method} \
			        --gene-lfc-method ${parameters.lfc_method} \
			        --normcounts-to-file \
              ${variance_params} \
			        ${rra_params} \
			        ${cnv_file} \
			        ${cnv_cellline}
	    fi
	    """
	else
	    """
      IFS=',' read -ra CTRL <<< "${parameters.control}"
      ctrl_string=""
      for i in "\${CTRL[@]}"; do
        if [ "\${ctrl_string}" == "" ]; then
          ctrl_string_buff="\${i}_ref_${parameters.name}"
        else
          ctrl_string_buff=",\${i}_ref_${parameters.name}"
        fi
        ctrl_string="\${ctrl_string}\${ctrl_string_buff}"
      done

	    prefilter_counts.R \
	        ${counts} \
	        \${ctrl_string} \
          ${parameters.treatment} \
	        0 > counts_filtered.txt

	    VERSION=\$(mageck -v 2>&1 >/dev/null)

	    if [ \$VERSION = "0.5.5" ]; then

		    mageck test \
		        --output-prefix ${parameters.name} \
		        --count-table counts_filtered.txt \
		        --control-id \${ctrl_string} \
		        --treatment-id ${parameters.treatment} \
		        --norm-method none \
		        --adjust-method ${parameters.fdr_method} \
		        --gene-lfc-method ${parameters.lfc_method} \
		        --normcounts-to-file
	    else
		    mageck test \
		        --output-prefix ${parameters.name} \
		        --count-table counts_filtered.txt \
		        --control-id \${ctrl_string} \
		        --treatment-id ${parameters.treatment} \
		        --norm-method none \
		        --adjust-method ${parameters.fdr_method} \
		        --gene-lfc-method ${parameters.lfc_method} \
		        --normcounts-to-file \
            ${variance_params} \
		        ${rra_params} \
		        ${cnv_file} \
		        ${cnv_cellline}
	    fi
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
