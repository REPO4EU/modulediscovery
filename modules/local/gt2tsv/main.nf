process GT2TSV {
    
    container 'docker.io/tiagopeixoto/graph-tool:release-2.58'

    input:
    path gt_file
    
    output:
    file("${gt_file.baseName}.nodes.tsv")

    script:
    """
    gt_to_tsv.py --input $gt_file  --output ${gt_file.baseName}.nodes.tsv
    """
}