include { SHORTEST_PATHS          } from '../../../modules/local/shortest_paths'
include { PHEN_TO_GENE            } from '../../../modules/local/phen_to_gene'
include { PROXIMITY               } from '../../../modules/local/proximity'

workflow GT_PROXIMITY {

    take:                                   // Workflow inputs
    network
    modules
    shortest_paths
    drug_to_target

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    SHORTEST_PATHS(network, shortest_paths)
    ch_versions = ch_versions.mix(SHORTEST_PATHS.out.versions)

    PHEN_TO_GENE(modules)
    ch_versions = ch_versions.mix(PHEN_TO_GENE.out.versions)

    PROXIMITY(network, SHORTEST_PATHS.out.sp, drug_to_target, PHEN_TO_GENE.out.phen_gene)
    ch_versions = ch_versions.mix(PROXIMITY.out.versions)

    emit:
    proxout   = PROXIMITY.out.proxout
    versions = ch_versions
}

