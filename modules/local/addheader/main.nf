process ADDHEADER {
    
    container 'docker.io/quirinmanz/gt2biopax:0.1.0'

    input:
    path file
    val header
    
     output:
    file("${file.baseName}.withHeader.txt")

    script:
    """
    addheader.py --input $file  --header $header   --output_file ${file.baseName}.withHeader.txt
    """
}