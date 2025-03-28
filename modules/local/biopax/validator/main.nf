process BIOPAX_VALIDATOR {
    label 'process_single'

    container "docker.io/quirinmanz/biopax-validator:6.0.0-SNAPSHOT"

    input:
    path biopax_files
    val validate_online

    output:
    path "biopax-validator-report.html", emit: validation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/sh

    java -version
    echo Running BioPAX Validator...

    if [ "${validate_online}" = true ]; then
        echo "Validating online"
        java -jar /biopax-validator/biopax-validator-client.jar . biopax-validator-report.html notstrict
    else
        echo "Validating offline"
        java --add-opens java.base/java.lang=ALL-UNNAMED --add-opens java.base/java.lang.reflect=ALL-UNNAMED -javaagent:/biopax-validator/lib/spring-instrument-5.3.29.jar -Xmx${task.memory.toGiga().toInteger()}g -Dfile.encoding=UTF-8 -Djava.security.egd=file:/dev/./urandom -jar /biopax-validator/biopax-validator.jar . --output=biopax-validator-report.html --profile=notstrict
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1)
        biopax-validator: \$(ls /biopax-validator/lib/biopax-validator*.jar)
    END_VERSIONS
    """
}
