//
// Prepares the input for DIAMOnD and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { ROBUST           } from '../../../modules/local/robust/main'

workflow GT_ROBUST {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "robust")                                  // Convert gt file to diamond specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    ROBUST(ch_seeds, GRAPHTOOLPARSER.out.network.collect())      // Run diamond on parsed network
    ch_versions = ch_versions.mix(ROBUST.out.versions.first())             // Collect versions


    emit:
    module   = ROBUST.out.module       // channel: [ module ]              emit the module discovered by DIAMOnD
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
