include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { SCA               } from '../../../modules/local/sca/main'

workflow GT_SCA {
    take:
    ch_seeds
    ch_network

    main:
    ch_versions = Channel.empty()
    GRAPHTOOLPARSER(ch_network, "sca")
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    ch_sca_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, meta, network]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "sca"
            [meta, seeds, network]
        }

    SCA(ch_sca_input)
    ch_versions = ch_versions.mix(SCA.out.versions.first())

    emit:
    module = SCA.out.module
    versions = ch_versions
}
