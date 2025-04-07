include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { SCA           } from '../../../modules/local/sca/main'
include { MODULEPARSER      } from '../../../modules/local/moduleparser/main'

workflow GT_SCA {
    take:
    ch_seeds
    ch_network

    main:
    ch_versions = Channel.empty()
    GRAPHTOOLPARSER(ch_network, "sca")
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    SCA(ch_seeds, GRAPHTOOLPARSER.out.network.collect())
    ch_versions = ch_versions.mix(SCA.out.versions.first())


    ch_module_seeds = SCA.out.module
        .join(ch_seeds, failOnMismatch: true, failOnDuplicate: true)
        .map{meta, module, seeds ->
            def dup = meta.clone()
            dup.id = meta.id + ".sca"
            [dup, module, seeds]
        }

        MODULEPARSER(ch_network, "sca", ch_module_seeds)
        ch_version = ch_versions.mix(MODULEPARSER.out.versions.first())


        emit:
        module = MODULEPARSER.out.network
        versions = ch_versions
}
