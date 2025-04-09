process SHORTEST_PATHS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path (network)

    output:
    tuple val(meta), path (network), path ("${meta.id}.shortest_paths.pkl"), emit: sp
    path "versions.yml", emit: versions

    script:
    """
    shortest_paths.py ${network} "${meta.id}.shortest_paths.pkl"

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            networkx: \$(python -c "import networkx; print(networkx.__version__)")
    END_VERSIONS
    """
}
