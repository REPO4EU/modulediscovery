process CHECKINPUT {
    
    container 'docker.io/quirinmanz/gt2biopax:0.1.0'

    input:
    path file
    

    script:
    """
    check_header.py --input $file  
    """
}