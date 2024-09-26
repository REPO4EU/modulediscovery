//
// Prepares the input for DOMINO and runs the tool
//

include { PREFIXLINES       } from '../../../modules/local/prefixlines/main'
include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DOMINO_SLICER     } from '../../../modules/local/domino/slicer/main'
include { DOMINO_DOMINO     } from '../../../modules/local/domino/domino/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_DOMINO {                        // Define the subworkflow, usually starts with the main input file format (.gt)
    take:                                   // Workflow inputs
    ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network                              // channel: [ val(meta[id,network_id]), path(network) ]


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    PREFIXLINES(ch_seeds, "entrez.")                                        // DOMINO interprets entrez ids as integers, so they are prefixed

    GRAPHTOOLPARSER(ch_network, "domino")                                   // Convert gt file to domino specific format, including prefixes
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    DOMINO_SLICER(GRAPHTOOLPARSER.out.network)                              // Run the DOMINO preprocessing step on the parsed networks
    ch_versions = ch_versions.mix(DOMINO_SLICER.out.versions)               // Collect versions


    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network), path(slices) ]
    ch_domino_input = PREFIXLINES.out
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, network]}, by: 0)
        .combine(DOMINO_SLICER.out.slices.map{meta, slices -> [meta.network_id, slices]}, by: 0)
        .map{network_id, meta, seeds, network, slices -> [meta, seeds, network, slices]}

    DOMINO_DOMINO(ch_domino_input)                                          // Run DOMINO on the preprocessed network
    ch_versions = ch_versions.mix(DOMINO_DOMINO.out.versions.first())       // Collect versions

    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module), path(seeds), path(network) ]
    ch_module_parser_input = DOMINO_DOMINO.out.modules                                // Extract the module
        .join(ch_seeds, failOnMismatch: true, failOnDuplicate: true)                  // Join with seed files
        .map{meta, module, seeds -> [meta.network_id, meta, module, seeds]}           // Combine with networks
        .combine(ch_network.map{meta, network -> [meta.network_id, network]}, by: 0)
        .map{network_id, meta, module, seeds, network ->                              // Adjust id
            def dup = meta.clone()
            dup.amim = "domino"
            dup.id = meta.id + "." + dup.amim
            dup.module_id = dup.id
            [ dup, module, seeds, network ]
        }

    MODULEPARSER(ch_module_parser_input, "domino")                    // Convert module from diamond specific format to gt file
    ch_versions = ch_versions.mix(MODULEPARSER.out.versions.first())


    emit:
    module   = MODULEPARSER.out.module  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
