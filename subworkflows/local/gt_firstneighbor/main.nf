//
// Prepares the input for FIRSTNEIGHBOR and runs the tool
//

include { FIRSTNEIGHBOR     } from '../../../modules/local/firstneighbor/main'
include { SPD     } from '../../../modules/local/spd/main'

workflow GT_FIRSTNEIGHBOR {
    take:                                   // Workflow inputs
    ch_seeds                                // File with seed genes
    ch_network                              // File with network in gt format
    spd_cutoff                              // SPD network cutoff


    main:

    ch_versions = Channel.empty()                                         // For collecting tool versions


    FIRSTNEIGHBOR(ch_seeds, ch_network)                                   // Run first neighbor
    ch_versions = ch_versions.mix(FIRSTNEIGHBOR.out.versions.first())     // Collect versions

    SPD(                                                                  // Run DOMINO
        FIRSTNEIGHBOR.out.module.collect(),                               // First input is the firstneighbor output module
        ch_network,                                                       // Second input is the input network
        spd_cutoff                                                        // Third input is the SPD cut-off
    )
    ch_versions = ch_versions.mix(SPD.out.versions.first())               // Collect versions


    emit:
    module_firstneighbor = FIRSTNEIGHBOR.out.module // channel: [ module ]  emit the module extracted using first neighbors
    module_spd = SPD.out.module         // channel: [ module ]              emit the module extracted using SPD refinement
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
