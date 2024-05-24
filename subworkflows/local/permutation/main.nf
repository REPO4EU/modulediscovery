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
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.original_seeds_id = meta.id
            [ dup, permuted_seeds]
        }
        .transpose()
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.id = permuted_seeds.baseName
            [ dup, permuted_seeds]
        }
    ch_versions = ch_versions.mix(SEEDPERMUTATION.out.versions)

    // Network expansion tools
    NETWORKEXPANSION(ch_permuted_seeds, ch_network)
    ch_permuted_modules = NETWORKEXPANSION.out.modules
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Evaluation
    // Input: initial module, initial seed genes, [randomized seeds, randomized modules]
    // Output:

    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
