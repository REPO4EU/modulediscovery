//
// Runs network permutation based evaluation of network expansion methods
//

include { NETWORKEXPANSION           } from '../networkexpansion'

workflow GT_NETWORKPERMUTATION {
    take:
    ch_modules  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]


    main:

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
