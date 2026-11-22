/*
 * Dorado is distributed as a standalone binary (not on Bioconda), so it must be on PATH
 * or provided through --dorado_bin. Basecalling with modified bases needs a GPU.
 */

process DORADO_BASECALLER {
    label 'process_gpu'

    input:
    path pod5_dir
    path reference

    output:
    path 'calls.bam', emit: bam

    script:
    // Offline/air-gapped runs: pass downloaded modification models as paths with --dorado_mod_models
    def mods = params.dorado_mod_models
        ? "--modified-bases-models ${params.dorado_mod_models}"
        : (params.dorado_modified_bases ? "--modified-bases ${params.dorado_modified_bases}" : '')
    """
    ${params.dorado_bin} basecaller \\
        ${params.dorado_model} \\
        ${pod5_dir} \\
        --recursive \\
        --device ${params.dorado_device} \\
        --kit-name ${params.dorado_kit} \\
        ${mods} \\
        --reference ${reference} \\
        ${params.dorado_args instanceof String ? params.dorado_args : ''} \\
        > calls.bam
    """

    stub:
    "touch calls.bam"
}

process DORADO_DEMUX {
    label 'process_medium'

    input:
    path calls_bam

    output:
    path 'barcodes/**.bam', emit: bams

    script:
    """
    ${params.dorado_bin} demux \\
        --no-classify \\
        --threads ${task.cpus} \\
        --output-dir barcodes \\
        ${calls_bam}
    """

    stub:
    """
    mkdir -p barcodes
    touch barcodes/stub_barcode01.bam barcodes/stub_barcode02.bam
    """
}
