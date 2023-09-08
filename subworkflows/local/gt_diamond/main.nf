//
// Check input samplesheet and get read channels
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { DIAMOND           } from '../../../modules/local/diamond/main'

workflow GT_DIAMOND {
    take:
    ch_seeds
    ch_network
    n                                      // DIAMOnD specific parameter "n"
    alpha                                  // DIAMOnD spefific parameter "alpha"


    main:

    ch_versions = Channel.empty()

    GRAPHTOOLPARSER(ch_network, "diamond")
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions.first())

    DIAMOND(ch_seeds, GRAPHTOOLPARSER.out.network, n, alpha)
    ch_versions = ch_versions.mix(DIAMOND.out.versions.first())


    emit:
    module   = DIAMOND.out.module                                     // channel: [ module ]
    versions = ch_versions                                            // channel: [ versions.yml ]
}
