process SPD {
    label 'process_single'

    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    path subnetwork
    path network
    val cutoff

    output:
    path "${subnetwork.baseName}_spd.gt", emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    spd_tool.py -s $subnetwork -n $network -o "${subnetwork.baseName}_spd.gt" -c $cutoff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
    END_VERSIONS
    """
}
