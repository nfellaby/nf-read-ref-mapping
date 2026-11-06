#!/usr/bin/env nextflow

include { MINIMAP2_INDEX    } from '../modules/local/minimap'
include { MINIMAP2_ALIGN    } from '../modules/local/minimap'
include { BAM_STATS         } from '../modules/local/bam_stats'
include { SUMMARISE_MAPPING } from '../modules/local/summarise_mapping'

workflow MAP_READS {
    take:
    ch_samples      // [ meta, [ fastq_1, (fastq_2) ], reference ]

    main:

        ch_versions = channel.empty()

        //
        // One minimap2 index per reference and preset, shared by every sample that uses it
        //
        ch_index_input = ch_samples
            .map { meta, _reads, reference ->
                [ [ id: meta.ref_id, path: reference.toString(), preset: meta.preset ], reference ]
            }
            .unique()

        MINIMAP2_INDEX ( ch_index_input )
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

        ch_align_input = ch_samples
            .map { meta, reads, reference ->
                [ [ id: meta.ref_id, path: reference.toString(), preset: meta.preset ], meta, reads ]
            }
            .combine(MINIMAP2_INDEX.out.index, by: 0)
            .map { _key, meta, reads, index -> [ meta, reads, index ] }

        //
        // Map, sort and index
        //
        MINIMAP2_ALIGN ( ch_align_input )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

        //
        // Per-sample mapping statistics
        //
        BAM_STATS ( MINIMAP2_ALIGN.out.bam )
        ch_versions = ch_versions.mix(BAM_STATS.out.versions.first())

        //
        // Cross-sample tables and HTML report
        //
        ch_manifest = BAM_STATS.out.stats
            .map { meta, _files ->
                "${meta.sample}\t${meta.platform}\t${meta.ref_id}\t${meta.single_end ? 'single' : 'paired'}\t${meta.id}"
            }
            .collectFile(name: 'mapping_manifest.tsv', newLine: true, sort: true)

        SUMMARISE_MAPPING (
            ch_manifest,
            BAM_STATS.out.stats.map { _meta, files -> files }.flatten().collect()
        )
        ch_versions = ch_versions.mix(SUMMARISE_MAPPING.out.versions)

    emit:
    bam      = MINIMAP2_ALIGN.out.bam
    summary  = SUMMARISE_MAPPING.out.summary
    html     = SUMMARISE_MAPPING.out.html
    versions = ch_versions
}
