process SHORTEST_PATHS {
    label 'process_single'

    container 'docker.io/kerstingjohannes/modulediscovery:1.0.0'

    input:
    path network
    path shortest_paths
    //path seed // still don't get it, no matter what are the seeds, the network's shortest paths should be the same ?

    output:
    path "${shortest_paths.name != 'NO_FILE' ? shortest_paths : 'shortest_paths.pkl'}", emit: sp
    path "versions.yml", emit: versions

    script:
    def output_sp = shortest_paths.name != 'NO_FILE' ? "$shortest_paths" : 'shortest_paths.pkl'
    // Check if shortest_paths does not exist
    if (!file(output_sp).exists()) {
        // Generate the shortest paths file
        """
        shortest_paths.py ${network} ${output_sp}

        cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                networkx: \$(python -c "import networkx; print(networkx.__version__)")
        END_VERSIONS
        """
    } else {
        // If the shortest paths file exists, do nothing
        """
        echo "Using provided shortest_paths: ${output_sp}"
        """
    }
}
