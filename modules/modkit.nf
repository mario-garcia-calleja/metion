process MODKIT_PILEUP {
    tag "${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/pileup", mode: params.publish_mode
    conda "bioconda::ont-modkit=0.4.3 bioconda::htslib=1.21"

    input:
    tuple val(sample), path(bam), path(bai)
    path reference

    output:
    tuple val(sample), path("${sample}.bed.gz"), path("${sample}.bed.gz.tbi"), emit: bedmethyl
    path "${sample}.pileup.log", emit: log

    script:
    // No --cpg by default: bacterial 6mA (e.g. Dam, GATC) and 5mC/4mC (e.g. Dcm, CCWGG)
    // are not in CpG context. Use --pileup_args "--motif GATC 1" to restrict to motifs.
    """
    modkit pileup ${bam} ${sample}.bed \\
        --ref ${reference} \\
        --threads ${task.cpus} \\
        --log-filepath ${sample}.pileup.log \\
        ${params.pileup_args instanceof String ? params.pileup_args : ''}
    bgzip -@ ${task.cpus} ${sample}.bed
    tabix -p bed ${sample}.bed.gz
    """

    stub:
    "touch ${sample}.bed.gz ${sample}.bed.gz.tbi ${sample}.pileup.log"
}

process MODKIT_DMR {
    tag "${comparison}"
    label 'process_medium'
    publishDir "${params.outdir}/dmr", mode: params.publish_mode
    conda "bioconda::ont-modkit=0.4.3"

    input:
    tuple val(comparison), path(ctrl_bed), path(ctrl_tbi), path(trt_bed), path(trt_tbi)
    path reference

    output:
    tuple val(comparison), path("${comparison}.dmr.bed"), emit: dmr
    path "${comparison}.dmr.log", emit: log

    script:
    // 21839 is the ChEBI code dorado uses for 4mC; it is counted together with 5mC at C sites
    def bases = params.dmr_bases.tokenize(',').collect { "--base ${it.trim()}" }.join(' ')
    """
    modkit dmr pair \\
        -a ${ctrl_bed} \\
        -b ${trt_bed} \\
        --ref ${reference} \\
        ${bases} \\
        --threads ${task.cpus} \\
        --log-filepath ${comparison}.dmr.log \\
        -o ${comparison}.dmr.bed
    """

    stub:
    "touch ${comparison}.dmr.bed ${comparison}.dmr.log"
}
