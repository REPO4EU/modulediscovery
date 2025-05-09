/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { INPUTCHECK               } from '../modules/local/inputcheck/main'
include { GRAPHTOOLPARSER          } from '../modules/local/graphtoolparser/main'
include { NETWORKANNOTATION        } from '../modules/local/networkannotation/main'
include { SAVEMODULES              } from '../modules/local/savemodules/main'
include { VISUALIZEMODULES         } from '../modules/local/visualizemodules/main'
include { GT2TSV as GT2TSV_Modules } from '../modules/local/gt2tsv/main'
include { GT2TSV as GT2TSV_Network } from '../modules/local/gt2tsv/main'
include { DIGEST                   } from '../modules/local/digest/main'
include { MODULEOVERLAP            } from '../modules/local/moduleoverlap/main'
include { DRUGPREDICTIONS          } from '../modules/local/drugpredictions/main'
include { TOPOLOGY                 } from '../modules/local/topology/main'
include { DRUGSTONEEXPORT          } from '../modules/local/drugstoneexport/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GT_BIOPAX             } from '../subworkflows/local/gt_biopax/main'
include { NETWORKEXPANSION      } from '../subworkflows/local/networkexpansion/main'
include { GT_SEEDPERMUTATION    } from '../subworkflows/local/gt_seedpermutation/main'
include { GT_NETWORKPERMUTATION } from '../subworkflows/local/gt_networkpermutation/main'
include { GT_PROXIMITY          } from '../subworkflows/local/gt_proximity/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//                                                //Evaluation
include { GPROFILER2_GOST        } from '../modules/nf-core/gprofiler2/gost/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_modulediscovery_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Channel guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ch_seeds: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    meta.id:
        - unique identifier for the specific seed file (combination of seeds_id and network_id)
        - used for process tags and output file names
    meta.seeds_id:
        - identifier for the input seed file (network_id not included)
        - used for naming modules (together with meta.id from ch_network)
        - used for merging ch_seeds and ch_modules (only in combination with meta.network_id!)
    meta.network_id:
        - used for merging ch_seeds and ch_network
        - used for merging ch_seeds and ch_modules (only in combination with meta.seeds_id!)
    seeds:
        - path to the seed file

ch_network: [ val(meta[id,network_id]), path(network) ]
    meta.id:
        - unique identifier for the specific network file (network_id)
        - used for process tags and output file names
        - used of naming modules (together with meta.seeds_id from ch_seeds)
    meta.network_id:
        - identifier for the input network file
        - used for merging ch_seeds and ch_network
    network:
        - path to the network file

ch_modules: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    meta.id:
        - unique identifier for the specific module file (combination of seeds_id, network_id, and amim)
        - used for process tags and output file names
    meta.module_id:
        - currently the same as meta.id
    meta.seeds_id:
        - identifier for the input seed file (network_id not included)
        - used for merging ch_modules and ch_seeds (only in combination with meta.network_id!)
    meta.network_id:
        - identifier for the input network file
        - used for merging ch_modules and ch_seeds (only in combination with meta.seeds_id!)
    meta.amim:
        - identifier for the network expansion method
    module:
        - path to the module file

*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MODULEDISCOVERY {


    take:
    ch_seeds                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network              // channel: [ val(meta[id,network_id]), path(network) ]
    ch_shortest_paths       // channel: [ val(meta[id,network_id]), path(shortest_paths) ]
    ch_permuted_networks    // channel: [ val(meta[id,network_id]), [path(permuted_networks)] ]

    main:

    // Params
    id_space = Channel.value(params.id_space)
    validate_online = Channel.value(params.validate_online)

    if(params.run_proximity){
        proximity_dt = file(params.drug_to_target, checkIfExists:true)
    }

    // Channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Run network parser for  networks, supported by graph-tool
    GRAPHTOOLPARSER(ch_network, 'gt')
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(GRAPHTOOLPARSER.out.multiqc)
    ch_network_gt = GRAPHTOOLPARSER.out.network


    // Check input
    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds), path(network) ]
    ch_seeds_network = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(ch_network_gt.map{meta, network -> [meta.network_id, network]}, by: 0)
        .map{key, meta, seeds, network -> [meta, seeds, network]}

    INPUTCHECK(ch_seeds_network)
    ch_seeds = INPUTCHECK.out.seeds
    INPUTCHECK.out.removed_seeds | view {meta, path -> log.warn("Removed seeds from $meta.id. Check multiqc report.") }
    ch_seeds_multiqc = INPUTCHECK.out.multiqc
        .map{ meta, path -> path }
        .collectFile(name: 'input_seeds_mqc.tsv', keepHeader: true)
    ch_multiqc_files = ch_multiqc_files.mix(ch_seeds_multiqc)

    // Add seeds modules to module channel
    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module)]
    ch_modules = INPUTCHECK.out.seeds_module
        .map{meta, path ->
            def dup = meta.clone()
            dup.amim = "no_tool"
            dup.id = meta.id + "." + dup.amim
            dup.module_id = dup.id
            [ dup, path ]
        }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        NETWORK EXPANSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Network expansion tools
    NETWORKEXPANSION(ch_seeds, ch_network_gt)
    ch_modules = ch_modules.mix(NETWORKEXPANSION.out.modules) // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module)]
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)


    // Annotate with network properties
    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module), path(network) ]
    ch_module_network = ch_modules
        .map{ meta, module -> [meta.network_id, meta, module]}
        .combine(ch_network_gt.map{meta, network -> [meta.network_id, network]}, by: 0)
        .map{newtork_id, meta, module, network -> [meta, module, network]}

    NETWORKANNOTATION(ch_module_network)
    ch_modules = NETWORKANNOTATION.out.module
    ch_versions = ch_versions.mix(NETWORKANNOTATION.out.versions)


    // Save modules
    SAVEMODULES(ch_modules)
    ch_versions = ch_versions.mix(SAVEMODULES.out.versions)

    // Visualize modules
    if(!params.skip_visualization){
        VISUALIZEMODULES(ch_modules, params.visualization_max_nodes)
        ch_versions = ch_versions.mix(VISUALIZEMODULES.out.versions)
    }

    // Drugstone export
    if(!params.skip_drugstone_export){
        DRUGSTONEEXPORT(ch_modules, id_space)
        ch_versions = ch_versions.mix(DRUGSTONEEXPORT.out.versions)
        ch_multiqc_files = ch_multiqc_files
            .mix(DRUGSTONEEXPORT.out.link.map{ meta, path -> path }.collectFile(name: 'drugstone_link_mqc.tsv', keepHeader: true))
    }

    // Annotation and BIOPAX conversion
    if(!params.skip_annotation){
        if( params.id_space != "symbol" & params.id_space != "ensembl" ){
            GT_BIOPAX(ch_modules, id_space, validate_online)
            ch_versions = ch_versions.mix(GT_BIOPAX.out.versions)
        } else {
            log.warn("Skipping annotation and BioPAX conversion (currently only uniprot or entrez IDs)")
        }

    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODULE EVALUATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if(!params.skip_evaluation){

        GT2TSV_Modules(ch_modules)
        GT2TSV_Network(ch_network_gt)

        // channel: [ val(meta), path(nodes) ]
        ch_nodes = GT2TSV_Modules.out

        // Module overlap
        ch_overlap_input = ch_nodes
            .multiMap { meta, nodes ->
                ids: meta.id
                nodes: nodes
            }
        MODULEOVERLAP(
            ch_overlap_input.ids.collect(),
            ch_overlap_input.nodes.collect()
        )
        ch_multiqc_files = ch_multiqc_files.mix(MODULEOVERLAP.out)

        // Topology evaluation
        TOPOLOGY(ch_modules)
        ch_versions = ch_versions.mix(TOPOLOGY.out.versions)
        ch_toplogy_multiqc = TOPOLOGY.out.multiqc
            .map{ meta, path -> path }
            .collectFile(name: 'topology_mqc.tsv', keepHeader: true)
        ch_multiqc_files = ch_multiqc_files.mix(ch_toplogy_multiqc)

        // Overrepresentation analysis
        if(!params.skip_gprofiler){

            ch_gprofiler_input = ch_nodes
                .map{ meta, path -> [meta.network_id, meta, path]}
                .combine(GT2TSV_Network.out.map{meta, path -> [meta.id, path]}, by: 0)
                .multiMap{key, meta, nodes, network ->
                    nodes: [meta, nodes]
                    network: network
                }

            GPROFILER2_GOST (
                ch_gprofiler_input.nodes,
                [],
                ch_gprofiler_input.network
            )
            ch_versions = ch_versions.mix(GPROFILER2_GOST.out.versions)
        }

        // Digest
        if(!params.skip_digest){

            ch_digest_input = ch_nodes
                .map{ meta, path -> [meta.network_id, meta, path]}
                .combine(ch_network_gt.map{meta, path -> [meta.id, path]}, by: 0)
                .multiMap{key, meta, nodes, network ->
                    nodes: [meta, nodes]
                    network: network
                }

            DIGEST (ch_digest_input.nodes, id_space, ch_digest_input.network, id_space)
            ch_versions = ch_versions.mix(DIGEST.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(
                DIGEST.out.multiqc
                .map{ meta, path -> path }
                .collectFile(name: 'digest_mqc.tsv', keepHeader: true)
            )
        }

        // Seed permutation based evaluation
        if(params.run_seed_permutation){
            GT_SEEDPERMUTATION(
                ch_modules.filter{ meta, path -> meta.amim != "no_tool" }, // Filter out no_tool modules
                ch_seeds,
                ch_network_gt
            )
            ch_versions = ch_versions.mix(GT_SEEDPERMUTATION.out.versions)
            ch_multiqc_files = ch_multiqc_files
                .mix(GT_SEEDPERMUTATION.out.multiqc_summary)
                .mix(GT_SEEDPERMUTATION.out.multiqc_jaccard)
        }

        // Network permutation based evaluation
        if(params.run_network_permutation){
            GT_NETWORKPERMUTATION(
                ch_modules.filter{ meta, path -> meta.amim != "no_tool" }, // Filter out no_tool modules
                ch_seeds,
                ch_network_gt,
                ch_permuted_networks
            )
            ch_versions = ch_versions.mix(GT_NETWORKPERMUTATION.out.versions)
            ch_multiqc_files = ch_multiqc_files
                .mix(GT_NETWORKPERMUTATION.out.multiqc_summary)
                .mix(GT_NETWORKPERMUTATION.out.multiqc_jaccard)
        }

    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DRUG PRIORITIZATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if(!params.skip_drug_predictions){
        def valid_algorithms = ['trustrank', 'closeness', 'degree']

        // Split the algorithms and check if they are valid
        ch_algorithms_drugs = Channel
            .of(params.drugstone_algorithms.split(','))
            .filter { algorithm ->
                if (!valid_algorithms.contains(algorithm)) {
                    throw new IllegalArgumentException("Invalid algorithm: $algorithm. Must be one of: ${valid_algorithms.join(', ')}")
                }
                return true
            }

        ch_drugstone_input = SAVEMODULES.out.nodes_tsv
            .combine(ch_algorithms_drugs)
            .multiMap { meta, module, algorithm ->
                module: [meta, module]
                algorithm: algorithm
            }
        includeIndirectDrugs = Channel.value(params.includeIndirectDrugs).map{it ? 1 : 0}
        includeNonApprovedDrugs = Channel.value(params.includeNonApprovedDrugs).map{it ? 1 : 0}
        DRUGPREDICTIONS(ch_drugstone_input.module, id_space, ch_drugstone_input.algorithm, includeIndirectDrugs, includeNonApprovedDrugs, params.result_size)
        ch_versions = ch_versions.mix(DRUGPREDICTIONS.out.versions)
    }

    // Drug prioritization - Proximity
    if(params.run_proximity){
        GT_PROXIMITY(ch_network_gt, SAVEMODULES.out.nodes_tsv, ch_shortest_paths, proximity_dt)
        ch_versions = ch_versions.mix(GT_PROXIMITY.out.versions)
    }


    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'modulediscovery_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
