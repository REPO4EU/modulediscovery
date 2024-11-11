//
// Prepares the input for FIRSTNEIGHBOR and runs the tool
//

include { FIRSTNEIGHBOR     } from '../../../modules/local/firstneighbor/main'

workflow GT_FIRSTNEIGHBOR {
    take:
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]

    main:

    ch_versions = Channel.empty()

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network) ]
    ch_seeds_network = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(ch_network.map{meta, network -> [meta.network_id, network]}, by: 0)
        .map{network_id, meta, seeds, network -> [meta, seeds, network]}

    FIRSTNEIGHBOR(ch_seeds_network)
    ch_versions = ch_versions.mix(FIRSTNEIGHBOR.out.versions.first())

    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_module = FIRSTNEIGHBOR.out.module
        .map{meta, path ->
            def dup = meta.clone()
            dup.amim = "firstneighbor"
            dup.id = meta.id + "." + dup.amim
            dup.module_id = dup.id
            [ dup, path ]
        }

    emit:
    module = ch_module      // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions  // channel: [ versions.yml ]
}
