process MAGECK {
    tag { parameters.name }

    publishDir path: "${params.outputDir}/${parameters.name}",
               mode: 'copy',
               overwrite: true,
               saveAs: { filename ->
                   if (filename.indexOf(".log") > 0) "$filename"
                   else if (filename.indexOf(".normalized.txt") > 0) "$filename"
                   else null
               }

    input:
    tuple val(parameters), path(counts), path(cnv), path(ctrl)

    output:
    tuple val("${parameters.name}"), path('*.sgrna_summary.txt'), path('*.gene_summary.txt'), emit: results
    path '*.log', emit: logs
    path '*.normalized.txt', emit: normalized

    script:
    def variance_params = parameters.variance_estimation ? "--variance-estimation-samples ${parameters.variance_estimation}" : ''
    def rra_params = params.min_rra_window > 0 ? "--additional-rra-parameters '-p ${params.min_rra_window}'" : ''
    def cnv_file = file(params.cnv).exists() & parameters.cnv_correction != '' ? "--cnv-norm ${cnv}" : ""
    def cnv_cellline = file(params.cnv).exists() & parameters.cnv_correction != '' ? "--cell-line ${parameters.cnv_correction}" : ""
    def control_sgRNAs = parameters.gene_level_statistics == 'control_sgRNAs' & file(params.control_sgRNAs).exists() ? "--control-sgrna ${ctrl}" : ''
    def estimate_min_count_from_samples = params.estimate_min_count_from_samples ? 1 : 0

    def filter = parameters.filter != "" & parameters.filter != null ? parameters.filter : parameters.control
    def control = parameters.control != "" ? parameters.control : "empty"
    def treatment = parameters.treatment != "" ? parameters.treatment : "empty"
    def variance = parameters.variance_estimation != "" & parameters.variance_estimation != null ? parameters.variance_estimation : "empty"

    if (parameters.norm_method == "quantile") {
        """
        prefilter_counts.R \
            ${counts} \
            ${filter} \
            ${params.min_count} \
            ${estimate_min_count_from_samples} \
            ${control} \
            ${treatment} \
            ${variance}  > counts_filtered.txt

        quantile_normalize_counts.R \
            counts_filtered.txt > counts_quantile_normalized.txt

        VERSION=\$(mageck -v 2>&1 >/dev/null)

        if [ \$VERSION = "0.5.5" ]; then
            mageck test \
                --output-prefix ${parameters.name} \
                --count-table counts_quantile_normalized.txt \
                --control-id ${parameters.control} \
                --treatment-id ${parameters.treatment} \
                --norm-method none \
                --adjust-method ${parameters.fdr_method} \
                --gene-lfc-method ${parameters.lfc_method} \
                --normcounts-to-file \
                ${control_sgRNAs}
        else
            mageck test \
                --output-prefix ${parameters.name} \
                --count-table counts_quantile_normalized.txt \
                --control-id ${parameters.control} \
                --treatment-id ${parameters.treatment} \
                --norm-method none \
                --adjust-method ${parameters.fdr_method} \
                --gene-lfc-method ${parameters.lfc_method} \
                --normcounts-to-file \
                ${control_sgRNAs} \
                ${variance_params} \
                ${rra_params} \
                ${cnv_file} \
                ${cnv_cellline}
        fi
        """
    } else {
        """
        prefilter_counts.R \
            ${counts} \
            ${filter} \
            ${params.min_count} \
            ${estimate_min_count_from_samples} \
            ${control} \
            ${treatment} \
            ${variance}  > counts_filtered.txt

        VERSION=\$(mageck -v 2>&1 >/dev/null)

        if [ \$VERSION = "0.5.5" ]; then
            mageck test \
                --output-prefix ${parameters.name} \
                --count-table counts_filtered.txt \
                --control-id ${parameters.control} \
                --treatment-id ${parameters.treatment} \
                --norm-method ${parameters.norm_method} \
                --adjust-method ${parameters.fdr_method} \
                --gene-lfc-method ${parameters.lfc_method} \
                --normcounts-to-file \
                ${control_sgRNAs}
        else
            mageck test \
                --output-prefix ${parameters.name} \
                --count-table counts_filtered.txt \
                --control-id ${parameters.control} \
                --treatment-id ${parameters.treatment} \
                --norm-method ${parameters.norm_method} \
                --adjust-method ${parameters.fdr_method} \
                --gene-lfc-method ${parameters.lfc_method} \
                --normcounts-to-file \
                ${control_sgRNAs} \
                ${variance_params} \
                ${rra_params} \
                ${cnv_file} \
                ${cnv_cellline}
        fi
        """
    }
}