#!/usr/bin/env nextflow
/*
 * map-reads: map Illumina and/or ONT reads to a per-sample reference with minimap2
 * and summarise how mapped reads are distributed across the reference.
 */
nextflow.enable.dsl = 2

include { MAP_READS } from './subworkflows/map_reads'

workflow {

    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }

    def input_ext = params.input.tokenize('.').last().toLowerCase()
    if (!(input_ext in ['csv', 'tsv'])) {
        error "Samplesheet must be .csv or .tsv: ${params.input}"
    }
    def sep = input_ext == 'tsv' ? '\t' : ','

    def presets = [ illumina: params.preset_illumina, ont: params.preset_ont ]

    //
    // Parse samplesheet: sample,platform,fastq_1,fastq_2,reference
    //
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true, sep: sep)
        .map { row ->
            if (!row.sample || !row.platform || !row.fastq_1 || !row.reference) {
                error "Each samplesheet row needs 'sample', 'platform', 'fastq_1' and 'reference' values: ${row}"
            }
            def platform = row.platform.toLowerCase()
            if (!(platform in presets.keySet())) {
                error "Unknown platform '${row.platform}' for ${row.sample}. Use illumina or ont."
            }
            def reference  = file(row.reference.trim(), checkIfExists: true)
            def single_end = !row.fastq_2?.trim()
            def meta = [
                id:         "${row.sample}_${platform}",
                sample:     row.sample,
                platform:   platform,
                single_end: single_end,
                preset:     presets[platform],
                ref_id:     reference.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
            ]
            def reads = single_end
                ? [ file(row.fastq_1.trim(), checkIfExists: true) ]
                : [ file(row.fastq_1.trim(), checkIfExists: true), file(row.fastq_2.trim(), checkIfExists: true) ]
            [ meta, reads, reference ]
        }

    // Each sample/platform pair must be unique, as it names the output files
    ch_samples
        .map { meta, _reads, _reference -> meta.id }
        .collect()
        .map { ids ->
            def dups = ids.countBy { id -> id }.findAll { _id, n -> n > 1 }.keySet()
            if (dups) {
                error "Duplicate sample/platform rows in samplesheet: ${dups.join(', ')}"
            }
        }

    MAP_READS ( ch_samples )

    MAP_READS.out.versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}
