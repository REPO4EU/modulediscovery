//
// Prepares the input for FIRSTNEIGHBOR and runs the tool
//

include { FIRSTNEIGHBOR     } from '../../../modules/local/firstneighbor/main'

workflow GT_FIRSTNEIGHBOR {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    GRAPHTOOLPARSER(ch_network, "firstneighbor")                            // Collect gt network file
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    FIRSTNEIGHBOR(ch_seeds, GRAPHTOOLPARSER.out.network.collect())          // Run first neighbor
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions.first())     // Collect versions


    emit:
    module   = FIRSTNEIGHBOR.out.module // channel: [ module ]              emit the module extracted using first neighbor
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
