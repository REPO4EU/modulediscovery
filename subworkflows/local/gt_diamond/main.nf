//
// Prepares the input for DIAMOnD and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DIAMOND           } from '../../../modules/local/diamond/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_DIAMOND {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format
    n                                       // DIAMOnD specific parameter "n"
    alpha                                   // DIAMOnD spefific parameter "alpha"


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "diamond")                                  // Convert gt file to diamond specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    // combine seed and network channel [meta, seeds, network]
    ch_seeds_network = ch_seeds
        .map{ meta, path -> [meta.network_id, meta, path]}
        .combine(GRAPHTOOLPARSER.out.network.map{meta, path -> [meta.id, path]}, by: 0)
        .map{key, meta, seeds, network -> [meta, seeds, network]}

    DIAMOND(ch_seeds_network, n, alpha)                                     // Run diamond on parsed network
    ch_versions = ch_versions.mix(DIAMOND.out.versions.first())

    ch_module_parser_input = DIAMOND.out.module                                    // Extract the module
        .join(ch_seeds_network, failOnMismatch: true, failOnDuplicate: true)        // Join with seed files
        .map{meta, module, seeds, network ->                                         // Adjust id
            def dup = meta.clone()
            dup.id = meta.id + ".diamond"
            dup.amim = "diamond"
            dup.seeds_id = meta.id
            [ dup, module, seeds, network ]
        }

    MODULEPARSER(ch_module_parser_input, "diamond")                    // Convert module from diamond specific format to gt file
    ch_versions = ch_versions.mix(MODULEPARSER.out.versions.first())


    emit:
    module   = MODULEPARSER.out.network       // channel: [ module ]              emit the module discovered by DIAMOnD
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
