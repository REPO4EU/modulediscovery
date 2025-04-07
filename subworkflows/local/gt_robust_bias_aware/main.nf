//
// Prepares the input for ROBUST_BIAS_AWARE and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { ROBUSTBIASAWARE   } from '../../../modules/local/robust_bias_aware/main'

workflow GT_ROBUSTBIASAWARE {
    take:
    ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network                              // channel: [ val(meta[id,network_id]), path(network) ]
    idspace

    main:

    ch_versions = Channel.empty()

    GRAPHTOOLPARSER(ch_network, "robust")
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    def idspaceUpper = idspace.toUpperCase()

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network) ]
    ch_robust_bias_aware_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, meta, network]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "robust_bias_aware"
            [meta, seeds, network]
        }

    ROBUSTBIASAWARE(ch_robust_bias_aware_input, idspaceUpper)
    ch_versions = ch_versions.mix(ROBUSTBIASAWARE.out.versions.first())



    emit:
    module   = ROBUSTBIASAWARE.out.module  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions                 // channel: [ versions.yml ]        emit collected versions
}
