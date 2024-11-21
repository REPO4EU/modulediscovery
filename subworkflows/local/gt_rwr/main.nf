//
// Prepares the input for RWR and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { RWR               } from '../../../modules/local/rwr/main'

workflow GT_RWR {
    take:                                   // Workflow inputs
    ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network                              // channel: [ val(meta[id,network_id]), path(network) ]
    scaling                                 // RWR specific parameter "scaling"
    symmetrical                             // RWR specific parameter "symmetrical"
    r                                       // RWR specific parameter "r"

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "rwr")                                      // Convert gt file to rwr specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network) ]
    ch_rwr_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, meta, network]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "rwr"
            [meta, seeds, network]
        }

    RWR(ch_rwr_input, scaling, symmetrical, r)                              // Run RWR on parsed network
    ch_versions = ch_versions.mix(RWR.out.versions.first())

    emit:
    module   = RWR.out.module  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
