//
// Runs network permutation based evaluation of network expansion methods
//

include { NETWORKEXPANSION             } from '../networkexpansion'
include { NETWORKPERMUTATION           } from '../../../modules/local/networkpermutation/main'
include { NETWORKPERMUTATIONEVALUATION } from '../../../modules/local/networkpermutationevaluation/main'

workflow GT_NETWORKPERMUTATION {
    take:
    ch_modules              // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    ch_seeds                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network              // channel: [ val(meta[id,network_id]), path(network) ]
    ch_permuted_networks    // channel: [ val(meta[id,network_id]), [path(permuted_networks)] ]

    main:
    ch_versions = Channel.empty()

    // Branch ch_permuted_networks based on whether the permuted_networks have already been pre-computed
    ch_permuted_networks = ch_network
        .join(ch_permuted_networks, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .branch{
            precomputed: it[2].size() > 0
            not_precomputed: true
        }

    // Permute the input network(s)
    // channel: [val(meta[id,network_id,n_permutations]), path(permuted_networks), val(output_name)]
    ch_permutation_input = ch_permuted_networks.not_precomputed
        // add n_permutations to meta
        .map{meta, network, permuted_networks -> [meta + [n_permutations: params.n_network_permutations], network]}
        // expand by number of permutations
        .combine(Channel.of(0..(params.n_network_permutations-1)))
        // add output file names
        .map{meta, network, permutation -> [meta, network, network.baseName + ".perm_" + permutation + "." + network.extension]}

    // Run network permutations
    NETWORKPERMUTATION(ch_permutation_input)
    ch_versions = ch_versions.mix(NETWORKPERMUTATION.out.versions)

    // Create required shape for NETWORKEXPANSION
    // channel: [val(meta[id,seeds_id,network_id,permuted_network_id,n_permutations]), path(permuted_network)]
    ch_permuted_networks =
        // Bring precomputed permuted networks in right shape
        ch_permuted_networks.precomputed.map{ meta, network, permuted_networks -> [meta, permuted_networks] }
        // Add n_permutations for later grouping
        .map{ meta, permuted_networks -> [ meta + [n_permutations: permuted_networks.size()], permuted_networks] }
        // Convert to long format
        .transpose()
        // Mix with precomputed permutations
        .mix(NETWORKPERMUTATION.out.permuted_network)
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
    ch_evaluation = ch_modules.map{meta, module -> [meta.module_id, meta, module]}
        // Join with permuted modules
        .join(
            ch_permuted_modules.map{meta, permuted_modules -> [meta.module_id, permuted_modules]},
            by: 0, failOnDuplicate: true, failOnMismatch: true
        )
        // Prepare channel for evaluation
        .multiMap{module_id, meta, module, permuted_modules ->
            module: [meta, module]
            permuted_modules: permuted_modules
        }

    NETWORKPERMUTATIONEVALUATION(
        ch_evaluation.module,
        ch_evaluation.permuted_modules,
    )
    ch_versions = ch_versions.mix(NETWORKPERMUTATIONEVALUATION.out.versions)
    ch_multiqc_summary =  NETWORKPERMUTATIONEVALUATION.out.multiqc_summary
        .map{ meta, path -> path }
        .collectFile(
            storeDir: "${params.outdir}/mqc_summaries",
            name: 'network_permutation_mqc.tsv',
            keepHeader: true
        )
    ch_multiqc_jaccard =
        NETWORKPERMUTATIONEVALUATION.out.multiqc_jaccard
        .map{ meta, path -> path }
        .collectFile(
            item -> "  " + item.text,
            storeDir: "${params.outdir}/mqc_summaries",
            name: 'network_permutation_jaccard_mqc.yaml',
            sort: true,
            seed: new File("$projectDir/assets/network_permutation_jaccard_header.yaml").text
        )



    emit:
    versions = ch_versions                  // channel: [ versions.yml ]        emit collected versions
    multiqc_summary = ch_multiqc_summary    // channel: [ multiqc_summary ]     emit collected multiqc files
    multiqc_jaccard = ch_multiqc_jaccard    // channel: [ multiqc_jaccard ]     emit collected multiqc files
}
