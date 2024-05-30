
include { SPD     } from '../../../modules/local/spd/main'

workflow GT_SPD {

    take:                                   // Workflow inputs
    ch_modules                              // Files with modules in gt format
    ch_network                              // File with network in gt format

    main:

    ch_versions = Channel.empty()                                         // For collecting tool versions

    SPD(
        ch_modules,                                                       // First input is the modules
        ch_network                                                        // Second input is the full network
    )
    ch_versions = ch_versions.mix(SPD.out.versions)                       // Collect versions

    ch_module_spd = SPD.out.module                                           // Extract the modules with SPD annotations
        .map{meta, path ->
            def dup = meta.clone()
            dup.id = meta.id + ".spd"
            [ dup, path ]
        }

    emit:
    module   = module_spd               // channel: [ module ]              emit the module containing SPD annotations
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
