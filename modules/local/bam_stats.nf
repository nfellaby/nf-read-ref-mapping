process BAM_STATS {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::samtools=1.25'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.{flagstat.json,coverage.tsv,mapq.tsv,windows.tsv}"), emit: stats
    path "versions.yml"                                                                , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools flagstat -O json $bam > ${prefix}.flagstat.json

    # Per-contig depth and breadth (samtools excludes unmapped, secondary, QC-fail and duplicates)
    samtools coverage $bam > ${prefix}.coverage.tsv

    # One pass over primary mapped reads (-F 0x904): MAPQ histogram and read starts per window
    touch ${prefix}.mapq.tsv ${prefix}.windows.tsv
    samtools view -F 0x904 $bam \\
        | awk -v w=${params.window_size} -v mq=${prefix}.mapq.tsv -v win=${prefix}.windows.tsv '
            { m[\$5]++; c[\$3 "\\t" int((\$4 - 1) / w)]++ }
            END {
                for (k in m) print k "\\t" m[k] > mq
                for (k in c) print k "\\t" c[k] > win
            }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
