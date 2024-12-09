
process NETWORKPERMUTATION {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(network), val(output_name)

    output:
    tuple val(meta), path("${output_name}"), emit: permuted_network
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    randomize_network.py --network ${network} --output ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
