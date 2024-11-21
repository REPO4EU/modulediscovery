//
// Prepares the input for DIAMOnD and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DIAMOND           } from '../../../modules/local/diamond/main'

workflow GT_DIAMOND {
    take:                                   // Workflow inputs
    ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network                              // channel: [ val(meta[id,network_id]), path(network) ]
    n                                       // DIAMOnD specific parameter "n"
    alpha                                   // DIAMOnD spefific parameter "alpha"


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "diamond")                                  // Convert gt file to diamond specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network) ]
    ch_diamond_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, meta, network]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "diamond"
            [meta, seeds, network]
        }

    DIAMOND(ch_diamond_input, n, alpha)                                     // Run diamond on parsed network
    ch_versions = ch_versions.mix(DIAMOND.out.versions.first())

    emit:
    module   = DIAMOND.out.module       // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
