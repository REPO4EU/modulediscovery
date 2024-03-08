process BIOPAX_PARSER {
    tag "$meta.id"
    label 'process_single'

//     conda "conda-forge::graph-tool=2.58"
    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    tuple val(meta), path(network)
    val idspace

    output:
    path "*.owl" , emit: biopax
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gt2biopax.py $network -i $idspace -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
        pybiopax: \$(python -c "import pybiopax; print(pybiopax.__version__)")
        nedrex: \$(python -c "import nedrex; print(nedrex.__version__)")
    END_VERSIONS
    """
}
