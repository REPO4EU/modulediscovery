include { GRAPHTOOLPARSER   } from '../../../modules/local/graphtoolparser/main'
include { TOPAS             } from '../../../modules/local/topas/main'

workflow GT_TOPAS {
    take:
    ch_seeds
    ch_network

    main:

    ch_version = Channel.empty()

    GRAPTHOOLPARSER(ch_network, "robust")
    ch_version = ch_version.mix(GRAPHTOOLPARSER.out.versions)

    ch_topas_input = ch_seeds
        .map{ meta, seeds -> [meta.network_id, meta, seeds]}
        .combine(GRAPHTOOLPARSER.out.network.map{ meta, network -> [meta.network_id, meta, network]}, by: 0)
        .map{network_id, seeds_meta, seeds, network_meta, network ->
            def meta = seeds_meta + network_meta
            meta.id = seeds_meta.seeds_id + "." + network_meta.id
            meta.amim = "topas"
            [meta, seeds, network]
        }
    TOPAS(ch_topas_input)
    ch_version = ch_version.mix(TOPAS.out.versions.first())

    emit:
    module   = TOPAS.out.module
    versions = ch_version
}