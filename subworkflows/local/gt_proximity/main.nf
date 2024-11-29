include { SHORTEST_PATHS          } from '../../../modules/local/shortest_paths'
include { PROXIMITY               } from '../../../modules/local/proximity'

workflow GT_PROXIMITY {

    take:                                   // Workflow inputs
    ch_network
    ch_modules
    shortest_paths
    drug_to_target

    main:

    ch_versions = Channel.empty()

    // Branch ch_shortest_paths based on whether the shortest paths have already been computed
    ch_shortest_paths = ch_network
        .map{meta, network ->
            [meta, network, file("${projectDir}/assets/NO_FILE", checkIfExists: true)]
        }
        .branch{
            no_sp: it[2].name == "NO_FILE"
            sp: true
        }

    // Compute shortest paths if they have not been computed
    // channel: [ val(meta[id,network_id]), path(network), path(sp) ]
    SHORTEST_PATHS(ch_shortest_paths.no_sp.map{meta, network, sp -> [meta, network]})
    ch_versions = ch_versions.mix(SHORTEST_PATHS.out.versions.first())

    // Combine the computed shortest paths with the input shortest paths
    // channel: [ val(meta[id,network_id]), path(network), path(sp) ]
    ch_shortest_paths = ch_shortest_paths.sp.mix(SHORTEST_PATHS.out.sp)

    //Prepare proximity input
    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module), path(network), path(sp)]
    ch_proximity_input = ch_modules
        .map{meta, module -> [meta.network_id, meta, module]}
        .combine(ch_shortest_paths.map{meta, network, sp -> [meta.network_id, network, sp]}, by: 0).view{it}
        .multiMap{network_id, meta, module, network, sp ->
            module: [meta, module]
            network: network
            sp: sp
        }

    PROXIMITY(
        ch_proximity_input.network,
        ch_proximity_input.sp,
        drug_to_target,
        ch_proximity_input.module
    )
    ch_versions = ch_versions.mix(PROXIMITY.out.versions.first())

    emit:
    //proxout   = PROXIMITY.out.proxout
    versions = ch_versions
}

