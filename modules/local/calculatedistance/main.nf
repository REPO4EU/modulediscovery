process CALCULATEDISTANCE {
    tag "$meta1.id"
    container "docker.io/quirinmanz/gt2biopax:0.1.1"

    input:
    tuple val(meta1), path(module)
    tuple val(meta2), path (nodes )

    output:
    tuple val(meta1), path("${meta1.id}.distance.txt")

    script:
    """
    calculate_distance.py --module $module --nodes $nodes --id ${meta1.id} --out ${meta1.id}.distance.txt
    """
}


