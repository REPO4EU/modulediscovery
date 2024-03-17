process DIGEST {

    container 'biocontainers/biodigest:0.2.16--pyhdfd78af_2'

    input:
    path target_file
    val target_type
    path network
    val network_type
    
    output:
  
    path("${target_file.baseName}")

    script:
    """
    digest.py --target_file $target_file  --target_type $target_type   --network $network  --network_type $network_type 
   
    """
}