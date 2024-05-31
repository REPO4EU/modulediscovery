//
// Runs input permutation based evluation of network expansion methods
//

include { NETWORKEXPANSION           } from '../networkexpansion'
include { SEEDPERMUTATION            } from '../../../modules/local/seedpermutation/main'

workflow PERMUTATION {
    take:
    ch_seeds                                // Files with seed genes
    ch_modules                              // Files with modules
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()

    // Permute the input seeds
    SEEDPERMUTATION(ch_seeds)
    ch_permuted_seeds = SEEDPERMUTATION.out.permuted_seeds
        // Add original meta.id as original_seeds_id
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.original_seeds_id = meta.id
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
    ch_versions = ch_versions.mix(SEEDPERMUTATION.out.versions)


    // Run network expansion tools on permuted seeds
    NETWORKEXPANSION(ch_permuted_seeds, ch_network)
    ch_permuted_modules = NETWORKEXPANSION.out.modules
        // Add original_seeds_id and amim to tuple for grouping
        .map{meta, modules ->
            [meta.original_seeds_id, meta.amim, modules]
        }
        // Group by original_seeds_id and amim
        .groupTuple(by: [0,1]).view{it[2].size()}
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Group network expansion ouputs for evaluation


    // Evaluation
    // Input: initial module, initial seed genes, [randomized seeds, randomized modules]
    // Output:

    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
