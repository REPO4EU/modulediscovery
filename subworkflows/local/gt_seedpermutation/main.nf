//
// Runs seed permutation based evaluation of network expansion methods
//

include { NETWORKEXPANSION           } from '../networkexpansion'
include { SEEDPERMUTATION            } from '../../../modules/local/seedpermutation/main'
include { SEEDPERMUTATIONEVALUATION      } from '../../../modules/local/seedpermutationevaluation/main'

workflow GT_SEEDPERMUTATION {
    take:
    ch_modules  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]


    main:

    ch_versions = Channel.empty()

    // Permute the input seeds
    SEEDPERMUTATION(ch_seeds)
    ch_versions = ch_versions.mix(SEEDPERMUTATION.out.versions)

    // Create required shape for NETWORKEXPANSION
    // channel: [val(meta[id,seeds_id,network_id,original_seeds_id,n_permutations]), path(permuted_seeds)]
    ch_permuted_seeds = SEEDPERMUTATION.out.permuted_seeds
        // Add original meta.id as original_seeds_id and n_permutations
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.original_seeds_id = meta.seeds_id
            dup.n_permutations = permuted_seeds.size()
            [ dup, permuted_seeds]
        }
        // Convert to long format
        .transpose()
        // Update id and seeds_id based on permuted seeds (original id is still stored as original_seeds_id)
        .map{meta, permuted_seeds ->
            def dup = meta.clone()
            dup.id = permuted_seeds.baseName + "." + meta.network_id
            dup.seeds_id = permuted_seeds.baseName
            [ dup, permuted_seeds]
        }


    // Run network expansion tools on permuted seeds
    NETWORKEXPANSION(ch_permuted_seeds, ch_network)
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Group by original_seeds_id, amim, and network_id to get one element per original module
    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), [path(permuted_modules)], [path(permuted_seeds)] ]
    ch_permuted_modules = NETWORKEXPANSION.out.modules
        //  Combine with permuted seeds
        .map{meta, permuted_module -> [meta.seeds_id, meta.network_id, meta, permuted_module]}
        .combine(ch_permuted_seeds.map{meta, permuted_seeds -> [meta.seeds_id, meta.network_id, permuted_seeds]}, by: [0,1])
        // Add original_seeds_id, amim, and network_id to tuple for grouping
        .map{seeds_id, network_id, meta, permuted_module, permuted_seeds ->
            key = groupKey(meta.subMap("original_seeds_id", "amim", "network_id"), meta.n_permutations)
            [key, meta, permuted_module, permuted_seeds]
        }
        // Group by original_seeds_id, amim, and network_id
        .groupTuple()
        // Add an ID (based on the original seeds)
        .map{key, meta, permuted_modules, permuted_seeds ->
            [ [ id: key.original_seeds_id + "." + key.network_id + "." + key.amim, module_id: key.original_seeds_id + "." + key.network_id + "." + key.amim, amim: key.amim, seeds_id: key.original_seeds_id, network_id: key.network_id], permuted_modules, permuted_seeds]
        }


    // Combine with original modules, seeds, and network
    // Shape: [val(meta[id,module_id,amim,seeds_id,network_id]), path(original_module), path(original_seeds), [path(permuted_modules)], [path(permuted_seeds)], network]
    ch_evaluation = ch_modules
        // Combine modules with seeds
        .map{meta, module -> [meta.seeds_id, meta.network_id, meta, module]}
        .combine(ch_seeds.map{meta, seeds -> [meta.seeds_id, meta.network_id, seeds]}, by: [0,1])
        .map{seeds_id, network_id, meta, module, seeds -> [meta.module_id, meta, module, seeds]}
        // Joine modules with permuted modules and seeds
        .join(ch_permuted_modules
            .map{ meta, permuted_modules, permuted_seeds ->
                [meta.module_id, permuted_modules, permuted_seeds]
            }, by: 0, failOnDuplicate: true, failOnMismatch: true
        )
        // Combine with network (key is network_id)
        .map{module_id, meta, module, seeds, permuted_modules, permuted_seeds ->
            [meta.network_id, meta, module, seeds, permuted_modules, permuted_seeds]
        }
        .combine(ch_network.map{meta, network-> [meta.network_id, network]}, by: 0)
        // Multimap to create the final shape
        .multiMap{network_id, meta, module, seeds, permuted_modules, permuted_seeds, network ->
            module: [meta, module]
            seeds: seeds
            permuted_seeds: permuted_seeds
            permuted_modules: permuted_modules
            network: network
        }


    // Evaluation
    SEEDPERMUTATIONEVALUATION(
        ch_evaluation.module,
        ch_evaluation.seeds,
        ch_evaluation.permuted_modules,
        ch_evaluation.permuted_seeds,
        ch_evaluation.network
    )
    ch_versions = ch_versions.mix(SEEDPERMUTATIONEVALUATION.out.versions)
    ch_multiqc_summary =  SEEDPERMUTATIONEVALUATION.out.multiqc_summary
        .map{ meta, path -> path }
        .collectFile(
            storeDir: "${params.outdir}/mqc_summaries",
            name: 'seed_permutation_mqc.tsv',
            keepHeader: true
        )
    ch_multiqc_jaccard =
        SEEDPERMUTATIONEVALUATION.out.multiqc_jaccard
        .map{ meta, path -> path }
        .collectFile(
            item -> "  " + item.text,
            storeDir: "${params.outdir}/mqc_summaries",
            name: 'seed_permutation_jaccard_mqc.yaml',
            sort: true,
            seed: new File("$projectDir/assets/seed_permutation_jaccard_header.yaml").text
        )



    emit:
    versions = ch_versions                  // channel: [ versions.yml ]        emit collected versions
    multiqc_summary = ch_multiqc_summary    // channel: [ multiqc_summary ]     emit collected multiqc files
    multiqc_jaccard = ch_multiqc_jaccard    // channel: [ multiqc_jaccard ]     emit collected multiqc files
}
