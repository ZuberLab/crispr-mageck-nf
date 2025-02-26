process POSTPROCESS {
    tag { name }

    publishDir path: "${params.outputDir}/${name}",
               mode: 'copy',
               overwrite: true

    input:
    tuple val(name), path(guides), path(genes), path(ctrl)

    output:
    tuple val(name), path('*_stats.txt'), emit: processed
    path '*.pdf', emit: qc

    script:
    def control_sgRNAs = file(params.control_sgRNAs).exists() ? "${ctrl}" : ''
    """
    postprocess_mageck.R ${guides} ${genes} ${control_sgRNAs}
    """
}