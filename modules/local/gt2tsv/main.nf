process GT2TSV {
    tag "$meta.id"
    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    tuple val(meta), path(gt_file)

    output:
    tuple val(meta), path("${meta.id}.nodes.tsv")

    script:
    """
    gt_to_tsv.py --input $gt_file  --output ${meta.id}.nodes.tsv
    """
}
