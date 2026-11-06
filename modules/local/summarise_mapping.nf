process SUMMARISE_MAPPING {
    label 'process_single'

    container 'community.wave.seqera.io/library/python:3.15.2--9307e9470900642f'
    
    input:
    path manifest
    path stats, stageAs: 'stats/*'

    output:
    path "mapping_summary.tsv"     , emit: summary
    path "contig_distribution.tsv" , emit: contigs
    path "window_distribution.tsv" , emit: windows
    path "mapq_distribution.tsv"   , emit: mapq
    path "mapping_summary.html"    , emit: html
    path "versions.yml"            , emit: versions

    script:
    """
    summarise_mapping.py \\
        --manifest $manifest \\
        --stats-dir stats \\
        --window-size ${params.window_size} \\
        --min-mapq ${params.min_mapq} \\
        --top-contigs ${params.top_contigs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
