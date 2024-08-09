process MODULEOVERLAP {
    label 'process_single'
    container 'docker.io/quirinmanz/gt2biopax:0.1.1'
    debug true

    input:
    val(ids)
    path(inputs)

    output:
    path('jaccard_similarity_matrix_mqc.tsv'), emit: jaccard_multiqc
    path('shared_nodes_matrix_mqc.tsv')      , emit: shared_multiqc

    script:
    """
    module_overlap.py --ids $ids --inputs $inputs
    """
}
