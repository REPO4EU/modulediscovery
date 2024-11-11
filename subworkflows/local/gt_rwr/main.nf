//
// Prepares the input for RWR and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { RWR               } from '../../../modules/local/rwr/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_RWR {
    take:                                   // Workflow inputs
    ch_seeds                                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network                              // channel: [ val(meta[id,network_id]), path(network) ]
    scaling                                 // RWR specific parameter "scaling"
    symmetrical                             // RWR specific parameter "symmetrical"
    r                                       // RWR specific parameter "r"

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "rwr")                                      // Convert gt file to rwr specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    // channel: [ val(meta[id,seeds_id,network_id), path(seeds), path(network) ]
    ch_rwr_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, network]}, by: 0)
        .map{network_id, meta, seeds, network -> [meta, seeds, network]}

    RWR(ch_rwr_input, scaling, symmetrical, r)                              // Run RWR on parsed network
    ch_versions = ch_versions.mix(RWR.out.versions.first())

    // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module), path(seeds), path(network) ]
    ch_module_parser_input = RWR.out.module                                           // Extract the module
        .join(ch_seeds, failOnMismatch: true, failOnDuplicate: true)                  // Join with seed files
        .map{meta, module, seeds -> [meta.network_id, meta, module, seeds]}           // Combine with networks
        .combine(ch_network.map{meta, network -> [meta.network_id, network]}, by: 0)
        .map{network_id, meta, module, seeds, network ->                              // Adjust id
            def dup = meta.clone()
            dup.amim = "rwr"
            dup.id = meta.id + "." + dup.amim
            dup.module_id = dup.id
            [ dup, module, seeds, network ]
        }

    MODULEPARSER(ch_module_parser_input, "rwr")                    // Convert module from rwr specific format to gt file
    ch_versions = ch_versions.mix(MODULEPARSER.out.versions.first())


    emit:
    module   = MODULEPARSER.out.module  // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
