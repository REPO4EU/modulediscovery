process GT2TSV {
    
    container "docker.io/quirinmanz/gt2biopax:0.1.0"

    input:
    path gt_file
    
    output:
    file("${gt_file.baseName}.nodes.tsv")

    script:
    """
    gt_to_tsv.py --input $gt_file  --output ${gt_file.baseName}.nodes.tsv
    """
}