//
// Prepares the input for FIRSTNEIGHBOR and runs the tool
//

include { FIRSTNEIGHBOR     } from '../../../modules/local/firstneighbor/main'

workflow GT_FIRSTNEIGHBOR {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format

    main:

    ch_versions = Channel.empty()                                         // For collecting tool versions

    // combine seed and network channel [meta, seeds, network]
    ch_seeds_network = ch_seeds
        .map{ meta, path -> [meta.network_id, meta, path]}
        .combine(ch_network.map{meta, path -> [meta.id, path]}, by: 0)
        .map{key, meta, seeds, network -> [meta, seeds, network]}

    FIRSTNEIGHBOR(ch_seeds_network)                                   // Run first neighbor
    ch_versions = ch_versions.mix(FIRSTNEIGHBOR.out.versions.first())     // Collect versions

    ch_module = FIRSTNEIGHBOR.out.module                       // Extract the module
        .map{meta, path ->
            def dup = meta.clone()
            dup.amim = "firstneighbor"
            dup.id = meta.id + "." + dup.amim
            [ dup, path ]
        }

    emit:
    module = ch_module // channel: [ module ]              emit the module extracted using first neighbors
    versions = ch_versions        // channel: [ versions.yml ]        emit collected versions
}
