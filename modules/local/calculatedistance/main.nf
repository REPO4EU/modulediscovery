process CALCULATEDISTANCE {
    tag "$meta.id"
    container "docker.io/quirinmanz/gt2biopax:0.1.1"

    input:
    tuple val(meta), path(module)

    output:
    tuple val(meta), path("${meta.id}.distance.tsv")

    script:
    """
    calculate_distance.py --module $module --id ${meta.id} --out ${meta.id}.distance.tsv
    """
}


