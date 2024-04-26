//
// Prepares the input for ROBUST_BIAS_AWARE and runs the tool
//

include { ROBUSTBIASAWARE } from '../../../modules/local/robust_bias_aware/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_ROBUSTBIASAWARE {
    take:
    ch_seeds
    ch_network
    idspace

    main:

    ch_versions = Channel.empty()

    def idspaceUpper = idspace.toUpperCase()
    ROBUSTBIASAWARE(ch_seeds, idspaceUpper)
    ch_versions = ch_versions.mix(ROBUSTBIASAWARE.out.versions.first())

    ch_module = ROBUSTBIASAWARE.out.module
        .map{meta, path ->
            def dup = meta.clone()
            dup.id = meta.id + ".robust_bias_aware"
            [ dup, path ]
        }

    MODULEPARSER(ch_network, "robust", ch_module)
    ch_versions = ch_versions.mix(MODULEPARSER.out.versions.first())


    emit:
    module   = MODULEPARSER.out.network
    versions = ch_versions
}
