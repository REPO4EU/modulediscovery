//
// Prepares the input for DOMINO and runs the tool
//

include { PREFIXLINES       } from '../../../modules/local/prefixlines/main'
include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DOMINO_SLICER     } from '../../../modules/local/domino/slicer/main'
include { DOMINO_DOMINO     } from '../../../modules/local/domino/domino/main'

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

    ch_domino_network = GRAPHTOOLPARSER.out.network
        .join(DOMINO_SLICER.out.slices, failOnMismatch: true, failOnDuplicate: true)

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network), path(slices) ]
    ch_domino_input = PREFIXLINES.out
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(ch_domino_network.map{ meta, network, slices -> [meta.network_id, meta, network, slices]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network, slices ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "domino"
            [meta, seeds, network, slices]
        }

    DOMINO_DOMINO(ch_domino_input)                                          // Run DOMINO on the preprocessed network
    ch_versions = ch_versions.mix(DOMINO_DOMINO.out.versions.first())       // Collect versions

    emit:
    module   = DOMINO_DOMINO.out.modules  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions                // channel: [ versions.yml ]        emit collected versions
}
