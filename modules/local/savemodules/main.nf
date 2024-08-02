process SAVEMODULES {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    tuple val(meta), path(module)

    output:
    tuple val(meta), path("${meta.id}.graphml"), emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    save_modules.py -m "${module}" -p "${meta.id}" -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
