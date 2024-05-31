//
// Prepares the input for RWR and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { RWR               } from '../../../modules/local/rwr/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_RWR {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format
    scaling                                 // RWR specific parameter "scaling"
    symmetrical                             // RWR specific parameter "symmetrical"
    r                                       // RWR specific parameter "r"

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "rwr")                                      // Convert gt file to rwr specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    RWR(ch_seeds, GRAPHTOOLPARSER.out.network.collect(), scaling, symmetrical, r)      // Run RWR on parsed network
    ch_versions = ch_versions.mix(RWR.out.versions.first())                 // Collect versions

    ch_module_seeds = RWR.out.module                                              // Extract the module
        .join(ch_seeds, failOnMismatch: true, failOnDuplicate: true)
        .map{meta, module, seeds ->
            def dup = meta.clone()
            dup.id = meta.id + ".rwr"
            dup.amim = "rwr"
            [ dup, module, seeds ]
        }

    MODULEPARSER(ch_network, "rwr", ch_module_seeds)


    emit:
    module   = MODULEPARSER.out.network
    versions = ch_versions
}
