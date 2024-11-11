//
// Runs network permutation based evaluation of network expansion methods
//

include { NETWORKEXPANSION           } from '../networkexpansion'
include { NETWORKPERMUTATION         } from '../../../modules/local/networkpermutation/main'

workflow GT_NETWORKPERMUTATION {
    take:
    ch_modules  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]

    main:
    ch_versions = Channel.empty()

    // Permute the input network(s)
    NETWORKPERMUTATION(ch_network)
    ch_versions = ch_versions.mix(NETWORKPERMUTATION.out.versions)

    // Create required shape for NETWORKEXPANSION
    // channel: [val(meta[id,seeds_id,network_id,original_network_id,n_permutations]), path(permuted_networks)]
    ch_permuted_networks = NETWORKPERMUTATION.out.permuted_networks
        // Add original meta.id as original_network_id and n_permutations
        .map{meta, permuted_networks ->
            def dup = meta.clone()
            dup.original_network_id = meta.network_id
            dup.n_permutations = permuted_networks.size()
            [ dup, permuted_networks]
        }
        // Convert to long format
        .transpose()
        // Update id and network_id based on permuted network (original id is still stored as original_network_id)
        .map{meta, permuted_network ->
            def dup = meta.clone()
            dup.id = permuted_network.baseName
            dup.network_id = dup.id
            [ dup, permuted_network]
        }
        //.view{it}

        // Generate module id based on seed_id and network_id
        // Network and seed id always belong to original? -> module id = meta.id + "." + meta.id -> return sum of maps?

    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
