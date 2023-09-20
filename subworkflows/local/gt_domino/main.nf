//
// Prepares the input for DOMINO and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DOMINO_SLICER     } from '../../../modules/local/domino/slicer/main'
include { DOMINO_DOMINO     } from '../../../modules/local/domino/domino/main'
include { PREFIXLINES       } from '../../../modules/local/prefixlines/main'

workflow GT_DOMINO {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    PREFIXLINES(ch_seeds, "entrez.")

    GRAPHTOOLPARSER(ch_network, "domino")                                   // Convert gt file to diamond specific format
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    DOMINO_SLICER(GRAPHTOOLPARSER.out.network.collect())                    // Run diamond on parsed network
    ch_versions = ch_versions.mix(DOMINO_SLICER.out.versions)               // Collect versions

    DOMINO_DOMINO(PREFIXLINES.out, GRAPHTOOLPARSER.out.network.collect(), DOMINO_SLICER.out.slices.collect())
    ch_versions = ch_versions.mix(DOMINO_DOMINO.out.versions.first())


    emit:
    module   = DOMINO_DOMINO.out.modules       // channel: [ modules ]             emit the module discovered by DIAMOnD
    versions = ch_versions                    // channel: [ versions.yml ]        emit collected versions
}
