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

    FIRSTNEIGHBOR(ch_seeds, ch_network)                                     // Run first neighbor
    ch_versions = ch_versions.mix(FIRSTNEIGHBOR.out.versions.first())       // Collect versions
    ch_module = FIRSTNEIGHBOR.out.module                                   // Extract the module
        .map{meta, path ->
            def dup = meta.clone()
            dup.id = meta.id + ".firstneighbor"
            [ dup, path ]
        }


    emit:
    module   = ch_module                // channel: [ module ]              emit the module extracted using first neighbor
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
