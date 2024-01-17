//
// Prepares the input for RWR and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { RWR               } from '../../../modules/local/rwr/main'

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


    emit:
    module   = RWR.out.module           // channel: [ module ]              emit the module discovered by RWR
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
