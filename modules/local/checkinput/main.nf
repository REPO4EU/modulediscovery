process CHECKINPUT {
    
    container 'docker.io/tiagopeixoto/graph-tool:release-2.58'

    input:
    path file
    

    script:
    """
    check_header.py --input $file  
    """
}