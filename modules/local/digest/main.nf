process DIGEST {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/biodigest:0.2.16--pyhdfd78af_2'

    input:
    tuple val(meta), path(target_file)
    val target_type
    path network
    val network_type

    output:
    tuple val(meta), path("${meta.id}"), emit: outdir
    path "versions.yml", emit: versions

    script:
    """
    digest.py --target_file $target_file  --target_type $target_type   --network $network  --network_type $network_type --outdir ${meta.id}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biodigest: \$(python -c "import digest; print(biodigest.__version__)")
    END_VERSIONS

    """
}
