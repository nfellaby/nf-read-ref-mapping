process MINIMAP2_INDEX {
    tag "${meta.id} (${meta.preset})"
    label 'process_high'

    conda 'bioconda::minimap2=2.31'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.mmi"), emit: index
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.${meta.preset}"
    """
    # The index stores the preset's k/w, so it is built once per reference and preset
    minimap2 \\
        -x ${meta.preset} \\
        -t ${task.cpus} \\
        $args \\
        -d ${prefix}.mmi \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::minimap2=2.31 bioconda::samtools=1.25'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    path "versions.yml"                                            , emit: versions

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def platform = meta.platform == 'ont' ? 'ONT' : 'ILLUMINA'
    def sort_threads = Math.max(1, task.cpus - 1)
    """
    minimap2 \\
        -a \\
        -x ${meta.preset} \\
        -t ${task.cpus} \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.sample}\\tPL:${platform}' \\
        $args \\
        $index \\
        $reads \\
        | samtools sort -@ ${sort_threads} -m 512M -o ${prefix}.bam -

    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}