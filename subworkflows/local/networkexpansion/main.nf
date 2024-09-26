//
// Runs all the expansion methods
//

include { GT_DIAMOND            } from '../gt_diamond'
include { GT_DOMINO             } from '../gt_domino'
include { GT_ROBUST             } from '../gt_robust'
include { GT_ROBUSTBIASAWARE    } from '../gt_robust_bias_aware'
include { GT_FIRSTNEIGHBOR      } from '../gt_firstneighbor'
include { GT_RWR                } from '../gt_rwr'

workflow NETWORKEXPANSION {
    take:
    ch_seeds    // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network  // channel: [ val(meta[id,network_id]), path(network) ]


    main:
    // Params
    diamond_n = Channel.value(params.diamond_n)
    diamond_alpha = Channel.value(params.diamond_alpha)

    rwr_scaling = Channel.value(params.rwr_scaling).map{it ? 1 : 0}
    rwr_symmetrical = Channel.value(params.rwr_symmetrical).map{it ? 1 : 0}
    rwr_r = Channel.value(params.rwr_r)

    id_space = Channel.value(params.id_space)


    ch_versions = Channel.empty()
    ch_modules  = Channel.empty()

    if(!params.skip_diamond){
        GT_DIAMOND(ch_seeds, ch_network, diamond_n, diamond_alpha)
        ch_versions = ch_versions.mix(GT_DIAMOND.out.versions)
        ch_modules = ch_modules.mix(GT_DIAMOND.out.module)
    }

    if(!params.skip_domino){
        GT_DOMINO(ch_seeds, ch_network)
        ch_versions = ch_versions.mix(GT_DOMINO.out.versions)
        ch_modules = ch_modules.mix(GT_DOMINO.out.module)
    }

    if(!params.skip_robust){
        GT_ROBUST(ch_seeds, ch_network)
        ch_versions = ch_versions.mix(GT_ROBUST.out.versions)
        ch_modules = ch_modules.mix(GT_ROBUST.out.module)
    }

    if(!params.skip_robust_bias_aware){
        GT_ROBUSTBIASAWARE(ch_seeds, ch_network, id_space)
        ch_versions = ch_versions.mix(GT_ROBUSTBIASAWARE.out.versions)
        ch_modules = ch_modules.mix(GT_ROBUSTBIASAWARE.out.module)
    }

    if(!params.skip_firstneighbor){
        GT_FIRSTNEIGHBOR(ch_seeds, ch_network)
        ch_versions = ch_versions.mix(GT_FIRSTNEIGHBOR.out.versions)
        ch_modules = ch_modules.mix(GT_FIRSTNEIGHBOR.out.module)
    }

    if(!params.skip_rwr){
        GT_RWR(ch_seeds, ch_network, rwr_scaling, rwr_symmetrical, rwr_r)
        ch_versions = ch_versions.mix(GT_RWR.out.versions)
        ch_modules = ch_modules.mix(GT_RWR.out.module)
    }


    emit:
    modules  = ch_modules               // channel: [ val(meta[id,module_id,amim,seeds_id,network_id]), path(module) ]
    versions = ch_versions              // channel: [ versions.yml ]
}
