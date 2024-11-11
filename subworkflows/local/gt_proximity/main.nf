include { SHORTEST_PATHS          } from '../../../modules/local/shortest_paths'
include { PHEN_TO_GENE            } from '../../../modules/local/phen_to_gene'
include { PROXIMITY               } from '../../../modules/local/proximity'

workflow GT_PROXIMITY {

    take:                                   // Workflow inputs
    ch_network
    ch_modules
    shortest_paths
    drug_to_target

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    SHORTEST_PATHS(ch_network, shortest_paths)
    ch_versions = ch_versions.mix(SHORTEST_PATHS.out.versions.first())

    PHEN_TO_GENE(ch_modules)
    ch_versions = ch_versions.mix(PHEN_TO_GENE.out.versions.first())

    ch_proximity_input = PHEN_TO_GENE.out.phen_gene
        .map{meta, module -> [meta.network_id, meta, module]}
        .combine(ch_network.map{meta, network -> [meta.network_id, network]}, by: 0)
        .combine(SHORTEST_PATHS.out.sp.map{meta, sp -> [meta.network_id, sp]}, by: 0)
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

