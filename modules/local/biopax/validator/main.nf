process BIOPAX_VALIDATOR {
    label 'process_single'

//     conda "conda-forge::graph-tool=2.58"
    container "docker.io/quirinmanz/biopax-validator:5.1.0"

    input:
    path biopax_files

    output:
    path "*.html", emit: validation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/sh

    java -version
    echo Running BioPAX Validator...

    java -javaagent:/biopax-validator/lib/spring-instrument-5.3.28.jar -Xmx2g -Dfile.encoding=UTF-8 --add-opens java.base/java.lang=ALL-UNNAMED -jar /biopax-validator/biopax-validator.jar . --output=biopax-validator-report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1)
        paxtools: \$(ls /biopax-validator/*javadoc.jar)
    END_VERSIONS
    """
}
