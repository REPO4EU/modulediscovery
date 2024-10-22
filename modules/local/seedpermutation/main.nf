
process SEEDPERMUTATION {
    tag "$meta.id"
    label 'process_single'
    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    tuple val(meta), path(seeds)

    output:
    tuple val(meta), path("${seeds.baseName}.*.${seeds.extension}") , emit: permuted_seeds
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    seed_permutation.py --seeds ${seeds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
