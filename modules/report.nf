process DMR_REPORT {
    tag "${comparison}"
    label 'process_low'
    publishDir "${params.outdir}/report/comparisons", mode: params.publish_mode
    conda "${projectDir}/environment.yml"

    input:
    tuple val(comparison), path(dmr)
    path regions

    output:
    path "${comparison}.region_summary.tsv", emit: summary
    path "${comparison}.stats.tsv",          emit: stats
    path "${comparison}.region_sites.tsv",   emit: sites
    path "${comparison}.significant_sites.tsv"

    script:
    def legacy = params.legacy ? '--legacy' : ''
    """
    metion_dmr.py \\
        --dmr ${dmr} \\
        --regions ${regions} \\
        --prefix ${comparison} \\
        --fdr ${params.fdr} \\
        --min-effect ${params.min_effect} \\
        --min-coverage ${params.min_coverage} \\
        --flank ${params.flank} \\
        ${legacy}
    """

    stub:
    """
    touch ${comparison}.region_summary.tsv ${comparison}.stats.tsv
    touch ${comparison}.region_sites.tsv ${comparison}.significant_sites.tsv
    """
}

process MERGE_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/report", mode: params.publish_mode
    conda "${projectDir}/environment.yml"

    input:
    path summaries
    path stats
    path sites

    output:
    path 'metion_summary.tsv'
    path 'metion_stats.tsv'
    path 'metion_report.xlsx'

    script:
    """
    metion_merge.py --summaries ${summaries} --stats ${stats} --sites ${sites}
    """

    stub:
    "touch metion_summary.tsv metion_stats.tsv metion_report.xlsx"
}
