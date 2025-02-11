process VISUALIZEMODULESDRUGS {
    tag "$meta_drugs.id"
    label 'process_single'

    input:
    tuple val(meta), path(module), val(meta_drugs), val(algorithm), path(drug_predictions)
    val max_nodes

    output:
    tuple val(meta_drugs), val(algorithm), path("${meta_drugs.id}.${algorithm}.pdf")  , emit: pdf  , optional: true
    tuple val(meta_drugs), val(algorithm), path("${meta_drugs.id}.${algorithm}.png")  , emit: png  , optional: true
    tuple val(meta_drugs), val(algorithm), path("${meta_drugs.id}.${algorithm}.svg")  , emit: svg  , optional: true
    tuple val(meta_drugs), val(algorithm), path("${meta_drugs.id}.${algorithm}.html") , emit: html , optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    visualize_modules.py -m "${module}" -p "${meta_drugs.id}.${algorithm}" -n ${max_nodes} -d ${drug_predictions} -l DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        graph-tool: \$(python -c "import graph_tool; print(graph_tool.__version__)")
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
        pyintergraph: \$(pip show pyintergraph | grep Version | awk '{print \$2}')
        pyvis: \$(pip show pyvis | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
