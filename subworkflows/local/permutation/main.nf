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

    // Create required shape for NETWORKEXPANSION
    // One channel element for each permuted seed file
    // Shape: [meta, permuted_seeds]
    ch_permuted_seeds = SEEDPERMUTATION.out.permuted_seeds
        // Add original meta.id as original_seeds_id and n_permutations
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.original_seeds_id = meta.id
            dup.n_permutations = permuted_seeds.size()
            [ dup, permuted_seeds]
        }
        // Convert to long format
        .transpose()
        // Update meta.id based on permuted seeds (original id is still stored as original_seeds_id)
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.id = permuted_seeds.baseName
            [ dup, permuted_seeds]
        }


    // Run network expansion tools on permuted seeds
    NETWORKEXPANSION(ch_permuted_seeds, ch_network)
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Group by original_seeds_id, amim, and network_id
    // One channel element for each original module
    // Shape: [meta, [permuted_modules], [permuted_seeds]]
    ch_permuted_modules = NETWORKEXPANSION.out.modules
        //  Combine with permuted seeds (key is seeds_id)
        .map{meta, permuted_module -> [meta.seeds_id, meta, permuted_module]}
        .combine(ch_permuted_seeds.map{meta, permuted_seeds -> [meta.id, permuted_seeds]}, by: 0)
        // Add original_seeds_id, amim, and network_id to tuple for grouping
        .map{seeds_id, meta, permuted_module, permuted_seeds ->
            key = groupKey(meta.subMap("original_seeds_id", "amim", "network_id"), meta.n_permutations)
            [key, permuted_module, permuted_seeds]
        }
        // Group by original_seeds_id, amim, and network_id
        .groupTuple()
        // Add an ID (based on the original seeds)
        .map{key, permuted_modules, permuted_seeds ->
            [ [ id: key.original_seeds_id + "." + key.amim, amim: key.amim, seeds_id: key.original_seeds_id, network_id: key.network_id], permuted_modules, permuted_seeds]
        }


    // Combine with original modules, seeds, and network
    // One channel element for each individual modules
    // Shape: [meta, original_module, original_seeds, [permuted_modules], [permuted_seeds], network]
    ch_evaluation = ch_modules
        // Combine modules with seeds (key is seeds_id)
        .map{meta, module -> [meta.seeds_id, meta, module]}
        .combine(ch_seeds.map{meta, seeds -> [meta.id, seeds]}, by: 0)
        .map{key, meta, module, seeds -> [meta, module, seeds]}
        // Combine with permuted modules and seeds (key is the entire meta map)
        .combine(ch_permuted_modules, by: 0)
        // Combine with network (key is network_id)
        .map{meta, module, seeds, permuted_modules, permuted_seeds ->
            [meta.network_id, meta, module, seeds, permuted_modules, permuted_seeds]
        }
        .combine(ch_network.map{meta, network-> [meta.id, network]}, by: 0)
        // Multimap to create the final shape
        .multiMap{network_id, meta, module, seeds, permuted_modules, permuted_seeds, network ->
            module: [meta, module]
            seeds: seeds
            permuted_seeds: permuted_seeds
            permuted_modules: permuted_modules
            network: network
        }


    // Evaluation
    PERMUTATIONEVALUATION(
        ch_evaluation.module,
        ch_evaluation.seeds,
        ch_evaluation.permuted_modules,
        ch_evaluation.permuted_seeds,
        ch_evaluation.network
    )
    ch_versions = ch_versions.mix(PERMUTATIONEVALUATION.out.versions)
    ch_multiqc_files = PERMUTATIONEVALUATION.out.multiqc_summary
        .map{ meta, path -> path }
        .collectFile(name: 'permutation_mqc.tsv', keepHeader: true)


    emit:
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
    multiqc_files = ch_multiqc_files          // channel: [ multiqc_summary ]     emit collected multiqc files
}
