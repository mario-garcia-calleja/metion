#!/usr/bin/env nextflow
/*
 * METION: bacterial DNA methylation from Oxford Nanopore reads
 *
 *   pod5 -> dorado basecaller (+mods, aligned) -> dorado demux -> modkit adjust/sort/index
 *        -> modkit pileup (bedMethyl) -> modkit dmr pair (per comparison)
 *        -> FDR filtering + region annotation -> merged report
 *
 * The pipeline can start from pod5 files, from per-sample BAMs or from bedMethyl files
 * (see the samplesheet columns in the README).
 */

include { DORADO_BASECALLER; DORADO_DEMUX } from './modules/dorado'
include { PREPARE_BAM                     } from './modules/samtools'
include { MODKIT_PILEUP; MODKIT_DMR       } from './modules/modkit'
include { DMR_REPORT; MERGE_REPORT        } from './modules/report'

def requireFile(String label, value) {
    if (!value) error "Missing required parameter --${label}"
    def f = file(value)
    if (!f.exists()) error "--${label}: file not found: ${value}"
    return f
}

// Relative paths in a samplesheet are resolved against the samplesheet's own folder
def sheetFile(sheet, String path) {
    def f = (path.startsWith('/') || path =~ /^[a-zA-Z0-9]+:\/\//) ? file(path) : file(sheet.parent.resolve(path).toString())
    if (!f.exists()) error "File listed in ${sheet.name} not found: ${path}"
    return f
}

workflow {
    log.info """
    M E T I O N   ~   bacterial methylation with Nanopore
    =====================================================
    samplesheet : ${params.samplesheet}
    comparisons : ${params.comparisons}
    pod5        : ${params.pod5 ?: '-'}
    reference   : ${params.reference}
    regions     : ${params.regions}
    outdir      : ${params.outdir}
    """.stripIndent()

    def reference   = requireFile('reference', params.reference)
    def regions     = requireFile('regions', params.regions)
    def samplesheet = requireFile('samplesheet', params.samplesheet)
    def comparisons = requireFile('comparisons', params.comparisons)

    samples = Channel.fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row without 'sample': ${row}"
            [
                sample   : row.sample,
                barcode  : row.barcode ? row.barcode.toString().padLeft(2, '0') : null,
                bam      : row.bam ? sheetFile(samplesheet, row.bam) : null,
                bedmethyl: row.bedmethyl ? sheetFile(samplesheet, row.bedmethyl) : null,
            ]
        }

    // ---- 1. bedMethyl supplied directly -------------------------------------------------
    from_bedmethyl = samples
        .filter { it.bedmethyl }
        .map { s -> tuple(s.sample, s.bedmethyl, file("${s.bedmethyl}.tbi", checkIfExists: true)) }

    // ---- 2. BAMs: supplied directly or produced from pod5 -------------------------------
    bams_given = samples
        .filter { !it.bedmethyl && it.bam }
        .map { s -> tuple(s.sample, s.bam) }

    if (params.pod5) {
        def by_barcode = samples
            .filter { !it.bedmethyl && !it.bam }
            .map { s ->
                if (!s.barcode) error "Sample ${s.sample}: needs a 'barcode' when starting from pod5"
                tuple(s.barcode, s.sample)
            }
        DORADO_BASECALLER(file(params.pod5, checkIfExists: true), reference)
        DORADO_DEMUX(DORADO_BASECALLER.out.bam)
        demuxed = DORADO_DEMUX.out.bams
            .flatten()
            .map { bam ->
                def m = bam.name =~ /barcode(\d+)/
                m ? tuple(m[0][1].padLeft(2, '0'), bam) : null
            }
            .filter { it != null }
            .groupTuple()                            // several BAM chunks per barcode
        bams_from_pod5 = by_barcode
            .join(demuxed, failOnMismatch: false)   // barcodes not in the samplesheet are dropped
            .map { barcode, sample, bams -> tuple(sample, bams) }
        all_bams = bams_given.mix(bams_from_pod5)
    } else {
        all_bams = bams_given
    }

    PREPARE_BAM(all_bams)
    MODKIT_PILEUP(PREPARE_BAM.out.bam, reference)
    bedmethyl = from_bedmethyl.mix(MODKIT_PILEUP.out.bedmethyl)

    // ---- 3. Pairwise comparisons -------------------------------------------------------
    pairs = Channel.fromPath(comparisons)
        .splitCsv(header: true, strip: true)
        .map { r -> tuple(r.control, r.comparison, r.treated) }
        .combine(bedmethyl, by: 0)                                   // add control files
        .map { ctrl, comp, trt, cbed, ctbi -> tuple(trt, comp, cbed, ctbi) }
        .combine(bedmethyl, by: 0)                                   // add treated files
        .map { trt, comp, cbed, ctbi, tbed, ttbi -> tuple(comp, cbed, ctbi, tbed, ttbi) }

    MODKIT_DMR(pairs, reference)
    DMR_REPORT(MODKIT_DMR.out.dmr, regions)
    MERGE_REPORT(
        DMR_REPORT.out.summary.collect(),
        DMR_REPORT.out.stats.collect(),
        DMR_REPORT.out.sites.collect(),
    )
}
