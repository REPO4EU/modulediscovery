process SHORTEST_PATHS {
    label 'process_single'

    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    path network
    path shortest_paths
    //path seed // still don't get it, no matter what are the seeds, the network's shortest paths should be the same ?

    output:
    path("*.{pkl,pickle,pcl}"), emit: sp
    path "versions.yml", emit: versions

    script:
    """
    shortest_paths.py ${network} ${shortest_paths}

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        networkx: \$(python -c "import networkx; print(networkx.__version__)")
END_VERSIONS
    """
}
