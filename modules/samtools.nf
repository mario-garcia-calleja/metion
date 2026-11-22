process PREPARE_BAM {
    tag "${sample}"
    label 'process_medium'
    conda "bioconda::samtools=1.21 bioconda::ont-modkit=0.4.3"

    input:
    tuple val(sample), path(bams, stageAs: 'input/*')

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam

    script:
    // Optionally drop modification codes before pileup (e.g. 'h' = 5hmC, as in the TFG)
    // An empty value on the command line arrives as boolean `true`: treat it, '' and 'none' as "keep all"
    def ignore = (params.ignore_mods instanceof String && params.ignore_mods && params.ignore_mods != 'none') ? params.ignore_mods : null
    def adjust = ignore
        ? "modkit adjust-mods combined.bam adjusted.bam --ignore ${ignore} --threads ${task.cpus}"
        : "ln -s combined.bam adjusted.bam"
    // dorado demux can write several BAMs per barcode: concatenate them first
    def combine = bams instanceof List && bams.size() > 1
        ? "samtools cat -o combined.bam ${bams.join(' ')}"
        : "ln -s ${bams instanceof List ? bams[0] : bams} combined.bam"
    """
    ${combine}
    ${adjust}
    samtools sort -@ ${task.cpus} -o ${sample}.bam adjusted.bam
    samtools index ${sample}.bam
    """

    stub:
    "touch ${sample}.bam ${sample}.bam.bai"
}
