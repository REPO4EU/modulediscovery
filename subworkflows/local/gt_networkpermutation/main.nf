//
// Runs network permutation based evaluation of network expansion methods
//

include { NETWORKEXPANSION             } from '../networkexpansion'
include { NETWORKPERMUTATION           } from '../../../modules/local/networkpermutation/main'
include { NETWORKPERMUTATIONEVALUATION } from '../../../modules/local/networkpermutationevaluation/main'

workflow GT_NETWORKPERMUTATION {
    take:
    ch_modules  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]

    main:
    ch_versions = Channel.empty()

    // Permute the input network(s)
    NETWORKPERMUTATION(ch_network, params.n_network_permutations)
    ch_versions = ch_versions.mix(NETWORKPERMUTATION.out.versions)

    // Create required shape for NETWORKEXPANSION
    // channel: [val(meta[id,seeds_id,network_id,permuted_network_id,n_permutations]), path(permuted_networks)]
    ch_permuted_networks = NETWORKPERMUTATION.out.permuted_networks
        // Add n_permutations
        .map{meta, permuted_networks ->
            def dup = meta.clone()
            dup.n_permutations = permuted_networks.size()
            [ dup, permuted_networks]
        }
        // Convert to long format
        .transpose()
        // Update id and permuted_network_id based on permuted network (original id is still stored as network_id) for the module parser
        .map{meta, permuted_network ->
            def dup = meta.clone()
            dup.id = permuted_network.baseName
            dup.permuted_network_id = dup.id
            [ dup, permuted_network]
        }

    // Run network expansion tools on permuted networks
    NETWORKEXPANSION(ch_seeds, ch_permuted_networks)
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Group by seeds_id, amim, and network_id to get one element per original module
    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), [path(permuted_modules)] ]
    ch_permuted_modules = NETWORKEXPANSION.out.modules
        .map{meta, permuted_module->
            key = groupKey(meta.subMap("seeds_id", "amim", "network_id"), meta.n_permutations)
            [key, meta, permuted_module]
        }
        // Group by seeds_id, amim, and network_id
        .groupTuple()
        // Add id and module_id based on the original modules
        .map{key, meta, permuted_modules ->
            [ [ id: key.seeds_id + "." + key.network_id + "." + key.amim, module_id: key.seeds_id + "." + key.network_id + "." + key.amim, amim: key.amim, seeds_id: key.seeds_id, network_id: key.network_id], permuted_modules]
        }

    // Join permuted modules with original modules
    ch_evaluation = ch_modules
        // Join with permuted modules
        .join(ch_permuted_modules, by: 0, failOnDuplicate: true, failOnMismatch: true)
        // Prepare channel for evaluation
        .multiMap{meta, module, permuted_modules ->
            module: [meta, module]
            permuted_modules: permuted_modules
        }

    NETWORKPERMUTATIONEVALUATION(
        ch_evaluation.module,
        ch_evaluation.permuted_modules,
    )
    ch_versions = ch_versions.mix(NETWORKPERMUTATIONEVALUATION.out.versions)


    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
