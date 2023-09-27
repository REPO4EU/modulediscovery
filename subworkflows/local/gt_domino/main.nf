//
// Prepares the input for DOMINO and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DOMINO_SLICER     } from '../../../modules/local/domino/slicer/main'
include { DOMINO_DOMINO     } from '../../../modules/local/domino/domino/main'
include { PREFIXLINES       } from '../../../modules/local/prefixlines/main'

workflow GT_DOMINO {                        // Define the subworkflow, usually starts with the main input file format (.gt)
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format


    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    PREFIXLINES(ch_seeds, "entrez.")                                        // DOMINO interprets entrez ids as integers, so they are prefixed

    GRAPHTOOLPARSER(ch_network, "domino")                                   // Convert gt file to domino specific format, including prefixes
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)             // Collect versions

    DOMINO_SLICER(GRAPHTOOLPARSER.out.network.collect())                    // Run the DOMINO preprocessing step on the parsed network
    ch_versions = ch_versions.mix(DOMINO_SLICER.out.versions)               // Collect versions

    DOMINO_DOMINO(                                                          // Run DOMINO
        PREFIXLINES.out,                                                    // First input is the prefixed seed gene sheet
        GRAPHTOOLPARSER.out.network.collect(),                              // Second input is the parsed network
        DOMINO_SLICER.out.slices.collect()                                  // Third input are the slices produced by the preprocessing step
    )
    ch_versions = ch_versions.mix(DOMINO_DOMINO.out.versions.first())       // Collect versions


    emit:                                                                   // Define, what the subworkflow will return
    module   = DOMINO_DOMINO.out.modules                                    // channel: [ modules ]             emit the modules discovered by DOMINO
    versions = ch_versions                                                  // channel: [ versions.yml ]        emit collected versions
}
