//
// Prepares the input for ROBUST and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { ROBUST           } from '../../../modules/local/robust/main'

workflow GT_ROBUST {
    take:
    ch_seeds
    ch_network


    main:

    ch_versions = Channel.empty()

    GRAPHTOOLPARSER(ch_network, "robust")
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    ROBUST(ch_seeds, GRAPHTOOLPARSER.out.network.collect())
    ch_versions = ch_versions.mix(ROBUST.out.versions.first())


    emit:
    module   = ROBUST.out.module
    versions = ch_versions
}
