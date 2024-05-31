//
// Runs input permutation based evluation of network expansion methods
//

include { NETWORKEXPANSION           } from '../networkexpansion'
include { SEEDPERMUTATION            } from '../../../modules/local/seedpermutation/main'
include { PERMUTATIONEVALUATION      } from '../../../modules/local/permutationevaluation/main'

workflow PERMUTATION {
    take:
    ch_seeds                                // Files with seed genes
    ch_modules                              // Files with modules
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()

    // Permute the input seeds
    SEEDPERMUTATION(ch_seeds)
    ch_versions = ch_versions.mix(SEEDPERMUTATION.out.versions)

    // Create shape [meta, permuted_seeds] for NETWORKEXPANSION
    ch_permuted_seeds = SEEDPERMUTATION.out.permuted_seeds
        // Add original meta.id as original_seeds_id and n_permutations
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.original_seeds = meta.id
            dup.n_permutations = permuted_seeds.size()
            [ dup, permuted_seeds]
        }
        // Convert to long format
        .transpose()
        // Update meta.id based on permuted seeds
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.id = permuted_seeds.baseName
            [ dup, permuted_seeds]
        }


    // Run network expansion tools on permuted seeds
    NETWORKEXPANSION(ch_permuted_seeds, ch_network)
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)


    // Create shape [meta, [permuted_modules], [permuted_seeds]]
    ch_permuted_modules = NETWORKEXPANSION.out.modules
        //  Combine with permuted seeds
        .map{meta, module -> [meta.seeds, meta, module]}
        .combine(ch_permuted_seeds.map{meta, seeds -> [meta.id, seeds]}, by: 0)
        // Add original_seeds and amim to tuple for grouping
        .map{seeds_id, meta, module, seeds ->
            key = groupKey(meta.subMap("original_seeds", "amim"), meta.n_permutations)
            [key, module, seeds]
        }
        // Group by original_seeds and amim
        .groupTuple()
        // Add an ID
        .map{key, modules, seeds ->
            [ [ id: key.original_seeds + "." + key.amim, amim: key.amim, seeds: key.original_seeds ], modules, seeds]
        }


    // Combine with original modules and seeds
    // Create shape [meta, original_module, original_seeds, [permuted_modules], [permuted_seeds]]
    ch_evaluation = ch_modules
        // Combine with original seeds
        .map{meta, module -> [[id: meta.seeds], meta, module]}
        .combine(ch_seeds, by: 0)
        .map{seeds_meta, meta, module, seeds -> [meta, module, seeds]}
        // Combine with original modules
        .combine(ch_permuted_modules, by: 0)
        .multiMap{meta, module, seeds, permuted_modules, permuted_seeds ->
            module: [meta, module]
            seeds: seeds
            permuted_seeds: permuted_seeds
            permuted_modules: permuted_modules
        }


    // Evaluation
    PERMUTATIONEVALUATION(
        ch_evaluation.module,
        ch_evaluation.seeds,
        ch_evaluation.permuted_modules,
        ch_evaluation.permuted_seeds
    )


    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
