//
// Prepares the input for ROBUST and runs the tool
//

include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { ROBUST            } from '../../../modules/local/robust/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

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

    ch_module_seeds = ROBUST.out.module
        .join(ch_seeds, failOnMismatch: true, failOnDuplicate: true)
        .map{meta, module, seeds ->
            def dup = meta.clone()
            dup.id = meta.id + ".robust"
            dup.amim = "robust"
            dup.seeds = meta.id
            [ dup, module, seeds ]
        }

    MODULEPARSER(ch_network, "robust", ch_module_seeds)
    ch_versions = ch_versions.mix(MODULEPARSER.out.versions.first())


    emit:
    module   = MODULEPARSER.out.network
    versions = ch_versions
}
